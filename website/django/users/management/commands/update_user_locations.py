from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from users.models import MyUser
from django.conf import settings
from argparse import FileType
import psycopg2
import csv
import requests
import pdb


# NOTE: To use this command, a Google Geocode API key must be set in site_settings.py or settings.py.
class Command(BaseCommand):
    help = 'Update user locations for existing users.'

    def handle(self, *args, **options):
        users = MyUser.objects.all()

        for user in users:
            address = getAddress(user.institution, user.city, user.state, user.country)

            GOOGLE_MAPS_API_URL = 'https://maps.googleapis.com/maps/api/geocode/json'

            try:
                params = {
                    'address': address,
                    'key': settings.GOOGLE_API_KEY
                }
            except NameError as err:
                print("Please add a GOOGLE_API_KEY variable with google geocode api key into settings.py before running.")
                sys.exit(1)

            # Do the request and get the response data
            req = requests.get(GOOGLE_MAPS_API_URL, params=params)
            res = req.json()

            # Use the first result
            try:
                result = res['results'][0]
                geodata = dict()
                user.latitude = result['geometry']['location']['lat']
                user.longitude = result['geometry']['location']['lng']
                user.save()
            except IndexError as err:
                print("User address undefined, setting lat/lng to empty strings.")
                user.latitude = ''
                user.longitude = ''
                user.save()


def addCommaIfNecessary(address):
    if len(address) > 0:
        address += ', '
    return address


def getAddress(institution, city, state, country):
    address = ''

    if len(city) > 0:
        address += city

    if len(state) > 0:
        address = addCommaIfNecessary(address)
        address += state

    if len(country) > 0:
        address = addCommaIfNecessary(address)
        address += country

    if len(address) < 4:
        address = institution

    if len(address) < 4:
        address = ""

    return address
