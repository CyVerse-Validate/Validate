# iPlant Agave API Python SDK

## 1. Introduction

This is the Python client SDK for the [iPlant Agave API](http://agaveapi.co). It is a pure Python library for interacting with Agave's RESTful services. All interaction with the Agave API requires valid authentication credentials. Agave uses [OAuth2](http://oauth.net/2) as its authorization mechanism. As a result, you will need a  valid bearer token to use this library. You can obtain a bearer token through a normal OAuth flow or you can obtain one directly from the [Agave API Store](https://agave.iplantc.org/store). 

## 2. Getting the code

The current version is **2.0.0-SNAPSHOT**. You can clone the library from git and include it directly in your project.

In order to run the unit tests, you will need to edit the `settings.xml` file by adding your client credentials. If you do not have a set of client credentials, please visit the [Agave API Store](https://agave.iplantc.org/store) to register a client application and get your credentials.

	> git clone https://bitbucket.org/taccaci/agave-sdk-python.git
	> cd agave-sdk-python


## 3. Agave Samples Demo

The following examples show how to programmatically run identical demos to those in the [Agave Samples](https://bitbucket.org/taccaci/agave-samples) project using the same Python client SDK.

### 3.0. Authentication

You can obtain an access token using your favorite OAuth2 library. If you are building an app that runs without a browser, then you can simply use your client and user credentials to obtain a token directly using the code below.

	import requests

	payload = {'grant_type': 'client_credentials', 'username': 'imicrobe', 'password': 'xxxxxxx', 'scope': 'PRODUCTION'}
	auth=('BOI3cMLqSweD21NxrwZbu3Eihmwa',' XXXXXXXXXXXXXXXX')
	r = requests.post('https://agave.iplantc.org/token', data=payload, auth=auth)
	print r.text

### 3.1. System Registration

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.2. File Management

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.3. App Listing

	require ("AgaveClient.python");

	$accessToken = 'abc123';
	$baseUrl = 'https://agave.iplantc.org';

	$apiClient = new AgaveClient($accessToken, $baseUrl);
	$appsService = AppsApi($apiClient);

	$response = $appsService->listApplication();

	foreach($response->result as $app) {
		print_r($app);
	}

### 3.4. Job Submission

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.5. Notification Management

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.6. Metadata Management

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.7. User Profile Discovery

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}

### 3.8. Internal User Management

	require ("AgaveClient.python");


	public class HelloWorld {

  		public static void main(String[] args) throws Exception {

  		}

  	}
