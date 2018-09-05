'''
Project: jdna
File: jbei_api
Author: Justin
Date: 2/9/17

Description: 

'''
import os
import json
import requests

class JBEI_API_Exception(Exception):
    """Generic Exception for JBEI API"""


class Verbose(object):
    def __init__(self):
        pass

    def __call__(self, f):
        def wrapped_f(obj, *args, **kwargs):
            print((f.__name__, 'started'))
            r = f(obj, *args, **kwargs)
            print((f.__name__, 'ended'))
            return r

        return wrapped_f


class RequestDecorator(object):
    def __init__(self, status_codes):
        if not isinstance(status_codes, list):
            status_codes = [status_codes]
        self.code = status_codes

    def __call__(self, f):
        def wrapped_f(*args, **kwargs):
            # Append home url
            args = list(args)
            args[1] = os.path.join(args[0].home, args[1])

            # Add default header
            default_header = {'Content-Type': 'application/json'}
            if not 'headers' in kwargs:
                kwargs['headers'] = {}
            default_header.update(kwargs['headers'])
            kwargs['headers'] = default_header

            # Make function call
            r = f(*args, **kwargs)

            # Raise error for status code
            if r.status_code not in self.code:
                http_codes = {
                    200: "OK - Request was successful",
                    201: "CREATED - Resource was created",
                    400: "BAD REQUEST",
                    401: "UNAUTHORIZED",
                    403: "FORBIDDEN",
                    404: "NOT FOUND",
                    500: "INTERNAL SERVER ERROR",
                    504: "SERVER TIMEOUT"}
                raise JBEI_API_Exception("HTTP Response Failed {} {}".format(
                    r.status_code, http_codes[r.status_code]))
            return r

        return wrapped_f


def JSON(f):
    def json_wrap(*args, **kwargs):
        return f(*args, **kwargs).json()

    return json_wrap


def Authenticate(f):
    def auth_wrap(*args, **kwargs):
        api = list(args)[0]

        if 'headers' not in kwargs:
            kwargs['headers'] = {}

        kwargs['headers'][JBEI_REST_API.SESSIONID_KEY] = api.sessionid
        r = f(*args, **kwargs)
        return r

    return auth_wrap


class JBEI_REST_API(object):
    SESSIONID_KEY = 'X-ICE-Authentication-SessionId'

    def __init__(self, email, password, home="https://public-registry.jbei.org/rest"):
        self.home = home
        self.auth = (email, password)
        self.accesstoken = self.create_access_token()
        self.sessionid = self.accesstoken['sessionId']

    @JSON
    def create_access_token(self):
        @RequestDecorator([200])
        def wrapped(self, x, **kwargs):
            return requests.post(x, **kwargs)

        return wrapped(self, 'accesstokens', data=json.dumps({'email': self.auth[0], 'password': self.auth[1]}))

    @RequestDecorator([200])
    @Authenticate
    def _post(self, what, **kwargs):
        pass


    @RequestDecorator([200])
    @Authenticate
    def _get(self, what, **kwargs):
        return requests.get(what, **kwargs)


    @RequestDecorator([200])
    @Authenticate
    def _delete(self, what, **kwargs):
        return requests.delete(what, **kwargs)


    @JSON
    def get_access_token(self):
        return self._get('accesstokens')


    def delete_access_token(self):
        return self._delete('accesstokens')


    @JSON
    def basic_search(self, query, **kwargs):
        #         default_search = dict(
        #             offset=0,
        #             limit=15,
        #             sort='relevance', # sort='RELEVANCE' | 'TYPE' | 'PART_ID' | 'CREATED',
        #             acs=True, # sort order
        #             searchWeb=False
        #         )
        search_params = {}
        search_params['query'] = query
        search_params.update(kwargs)
        return self._get('search', params=search_params)


    @JSON
    def advanced_search(self, searchWeb=False, **kwargs):
        default_search = {
            'queryString': '',  # String,
            'blastQuery': {
                'blastProgram': 'BLAST_N',  # 'BLAST_N'/'TBLAST_X',
                'sequence': ''  # String
            },
            'bioSafetyOption': '',  # 'LEVEL_ONE'/'LEVEL_TWO',
            'entryTypes': ['STRAIN', 'PLASMID', 'PART', 'ARABIDOPSIS'],
            'parameters': {
                'sortField': None,  # 'RELEVANCE', # 'RELEVANCE'/'TYPE'/'PART_ID'/'CREATED'/'NAME',
                'sortAscending': True,  # Boolean,
                'start': 0,  # Integer,
                'retrieveCount': 1000,  # Integer,
                'hasSequence': None,  # true/false/undefined,
                'hasAttachment': None,  # true/false/undefined,
                'hasSample': None,  # true/false/undefined
            },
            'fieldFilters': None  # [{'field': 'NAME'/'ALIAS','PART_ID', 'filter': String}]
        }

        def _clean_dictionary(dic):
            keys = list(dic.keys())
            for k in keys:
                v = dic[k]
                if v is None or v == '':
                    dic.pop(k)
            return dic

        default_search.update(kwargs)
        _clean_dictionary(default_search)
        if default_search['blastQuery']['sequence'] == '':
            default_search.pop('blastQuery')

        _clean_dictionary(default_search['parameters'])
        print(default_search)
        r = self._post('search', params={'searchWeb': searchWeb}, data=json.dumps(default_search))
        print((r.url))
        return r


jbei = JBEI_REST_API("email", "password")
# jbei.get_access_token()


print((len(jbei.basic_search('pBAD')['results'])))
# print jbei.advanced_search({'queryString': 'pBAD'})
# jbei.delete_access_token()

# search_params = {'blastQuery': {
#     'blastProgram': 'BLAST_N',  # 'BLAST_N'/'TBLAST_X',
#     'sequence': '',  # String}