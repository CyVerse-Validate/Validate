#!/usr/bin/env python
"""
Copyright 2012 Wordnik, Inc.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
"""
class MetadataSchema:
    """NOTE: This class is auto generated by the swagger code generator program.
    Do not edit the class manually."""


    def __init__(self):
        self.swaggerTypes = {
            'created': 'date-time',
            'internalUsername': 'str',
            'lastUpdated': 'date-time',
            'owner': 'str',
            'schema': 'str',
            'uuid': 'str'

        }


        #A timestamp indicating when this Metadata was created in the metadata schema store.
        self.created = None # date-time
        #The name of the Internal User, if any, who owns this schema.
        self.internalUsername = None # str
        #A timestamp indicating when this Metadata was last updated in the metadata schema store.
        self.lastUpdated = None # date-time
        #The API user who owns this Schema.
        self.owner = None # str
        #A JSON Schema
        self.schema = None # str
        #The UUID for this Schema.
        self.uuid = None # str
        