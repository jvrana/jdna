'''
Project: jdna
File: json_dna
Author: Justin
Date: 2/19/17

Description: 

'''

var generalizedJsonFormat = {
    "size" : 25,
    "sequence" : "asaasdgasdgasdgasdgasgdasgdasdgasdgasgdagasdgasdfasdfdfasdfa",
    "circular" : true,
    "name" : "pBbS8c-RFP",
    "description" : "",
    "features" : [
        {
            "name" : "anonymous feature",
            "type" : "misc_feature",
            "id" : "5590c1978979df000a4f02c7", //Must be a unique id. If no id is provided, we'll autogenerate one for you
            "start" : 1,
            "end" : 3,
            "strand" : 1,
            "notes" : {},
        },
        {
            "name" : "coding region 1",
            "type" : "CDS",
            "id" : "5590c1d88979df000a4f02f5",
            "start" : 12,
            "end" : 9,
            "strand" : -1,
            "notes" : {},
        }
    ],
}



{
  "aliases": [
    "pDA001"
  ],
  "annotations": [
    {
      "color": "#303030",
      "end": 25,
      "name": "Untitled",
      "start": 16,
      "strand": 1,
      "type": "gene"
    }
  ],
  "bases": "GATTACATA",
  "circular": true,
  "color": "#303030",
  "createdAt": "2014-06-30T19:52:33.798355+00:00",
  "creator": {
    "avatarUrl": "test/a/90YumRpRp5A0aD/ent_Debas123ldk-lpasteur.png",
    "email": "lpasteur@benchling.test",
    "handle": "lpasteur",
    "id": "ent_Debas123ldk",
    "name": "Louis Pasteur"
  },
  "description": "A short description",
  "editURL": "/test/f/tur36kHe-folder_one/seq-nIlhJqUW-sequence_one/edit",
  "folder": {
    "id": "lib_ABCdlekf12",
    "name": "Untitled"
  },
  "id": "seq_1243bka",
  "length": 9,
  "modifiedAt": "2014-06-30T19:52:33.798355+00:00",
  "name": "Untitled",
  "notes": [
    {
      "created_at": "2014-06-30T19:52:33.798355+00:00",
      "creator": "ent_Debas123ldk",
      "end": 25,
      "start": 16,
      "text": "Some text"
    }
  ],
  "primers": [
    {
      "bases": "GATTACATA",
      "bind_position": 16,
      "color": "#303030",
      "created_at": "2014-06-30T19:52:33.798355+00:00",
      "end": 25,
      "name": "Untitled",
      "overhang_length": 3,
      "start": 16,
      "strand": 1
    }
  ],
  "tagSchema": {
    "folderItemType": "sequence",
    "id": "ts_s2fjAlb13k",
    "name": "a tag schema"
  },
  "tags": [
    {
      "name": "Untitled",
      "reference": "seq_1243bka",
      "url": "/url/to/resource",
      "value": "a value"
    }
  ]
}