#!/usr/bin/env python                                                           

import argparse
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new", help="New(er) JSON file")
    parser.add_argument("--old", help="Old(er) JSON file")
    args = parser.parse_args()

    newJson = json.loads(open(args.new).read())
    oldJson = json.loads(open(args.old).read())
    addedEntries = dict()
    removedEntries = dict()
    changedEntries = dict()

    print "Present in new but not old"
    for  newKey in newJson.keys():
        if not oldJson.has_key(newKey):
            print newKey, newJson[newKey]
            newString = str(newJson[newKey])
            if not addedEntries.has_key(newString):
                addedEntries[newString] = list()
            addedEntries[newString].append(newKey)
    print "Added entries"
    for key in addedEntries.keys():
        print key, addedEntries[key]
    print "Present in old but not new"
    for  oldKey in oldJson.keys():
        if not newJson.has_key(oldKey):
            print oldKey, oldJson[oldKey]
            oldString = str(oldJson[oldKey])
            if not removedEntries.has_key(oldString):
                removedEntries[oldString] = list()
            removedEntries[oldString].append(oldKey)
    print "Removed entries"
    for key in removedEntries.keys():
        print key, removedEntries[key]
    print "Changed between old and new"
    for newKey in newJson.keys():
        if oldJson.has_key(newKey):
            if newJson[newKey] != oldJson[newKey]:
                print newKey
                print "\told:", oldJson[newKey]
                print "\tnew:", newJson[newKey]
            

if __name__ == "__main__":
    main()
            
