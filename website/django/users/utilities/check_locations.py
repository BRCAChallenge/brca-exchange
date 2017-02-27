#!/usr/bin/env python

import argparse
import csv
import requests
import pdb


API_KEY = "AIzaSyC2Poc7PB1X6TgJ06XdXRkSVagB_pN5gzw"
OUTPUT_FIELDNAMES = {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--users", default="users.csv",
                        help="CSV dump from users_myuser table in DB")
    parser.add_argument("--output", default="bad_locations.csv",
                        help="File with users that have bad locations.")

    args = parser.parse_args()
    users = csv.DictReader(open(args.users, "r"), delimiter=",")

    fieldnames = users.fieldnames

    fieldnames.extend(['g_address', 'g_latitude', 'g_longitude', 'correct_location'])

    with open(args.output, 'wb') as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        for user in users:
            user_id = user['id']
            institution = user['institution']
            city = user['city']
            state = user['state']
            country = user['country']
            lat = user['latitude']
            lng = user['longitude']

            # current address method
            address = "" + institution + "," + city + "," + state + "," + country

            # todo: test better address methods

            GOOGLE_MAPS_API_URL = 'https://maps.googleapis.com/maps/api/geocode/json'

            params = {
                'address': address,
                'key': API_KEY
            }

            # Do the request and get the response data
            req = requests.get(GOOGLE_MAPS_API_URL, params=params)
            res = req.json()

            # Use the first result
            try:
                result = res['results'][0]
                geodata = dict()
                geodata['lat'] = result['geometry']['location']['lat']
                geodata['lng'] = result['geometry']['location']['lng']
                geodata['address'] = result['formatted_address']

                user['correct_location'] = within_one(geodata['lat'], lat) and within_one(geodata['lng'], lng)

                user['g_address'] = geodata['address']
                user['g_latitude'] = geodata['lat']
                user['g_longitude'] = geodata['lng']

                try:
                    user = convert(user)
                    writer.writerow(user)
                except UnicodeEncodeError as err:
                    print "Error encoding"
                    pdb.set_trace()

            except IndexError as err:
                user['correct_location'] = user['latitude'] == '' and user['longitude'] == ''
                user['g_address'] = ''
                user['g_latitude'] = ''
                user['g_longitude'] = ''
                try:
                    user = convert(user)
                    writer.writerow(user)
                except UnicodeEncodeError as err:
                    print "Error encoding"
                    pdb.set_trace()


def within_one(first, second):
    try:
        return abs(float(first) - float(second)) <= 1
    except ValueError as err:
        if first is None:
            first = ''
        if second is None:
            second = ''
        if first != second:
            return False
        else:
            return True


def convert(input):
    if isinstance(input, dict):
        return {convert(key): convert(value) for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [convert(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('ascii', 'ignore')
    else:
        return input


if __name__ == "__main__":
    main()
