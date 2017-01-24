from httmock import all_requests, response, urlmatch, HTTMock
import pytest
import unittest
import requests
import os
import datetime
import json
import django
from brca import settings
import pdb

os.environ['DJANGO_SETTINGS_MODULE'] = 'brca.settings'
django.setup()

from django.test import TestCase
from django.utils import timezone
from django.test.client import RequestFactory
from django.contrib.auth.models import AnonymousUser
from .models import MyUser
from users import views

settings.MAILCHIMP_URL = 'https://mailchimp.com'
settings.MAILCHIMP_KEY = '12345'
settings.MAILCHIMP_LIST = '12345'


@urlmatch(netloc=r'(.*\.)?google\.com$')
def google_mock(url, request):
    headers = {'Set-Cookie': 'foo=bar;'}
    content = {'success': True}
    return response(200, content, headers, None, 5, request)


@all_requests
def response_content(url, request):
    headers = {'Set-Cookie': 'foo=bar;'}
    content = {'success': True}
    return response(200, content, headers, None, 5, request)


class TestViewsAPI(TestCase):
    def setUp(self):
        self.factory = RequestFactory()
        self.test_user = self.create_test_user()

    def test_mailinglist_signup_success(self):
        with HTTMock(response_content):
            request = self.factory.post('/mailinglist/', {
                                        'email': 'test@example.com',
                                        'firstName': 'test',
                                        'lastName': 'example'
                                        })
            response = views.mailinglist(request)
            self.assertEqual(response.status_code, 200)

    def test_mailinglist_signup_failure(self):
        with HTTMock(google_mock):
            request = self.factory.post('/mailinglist/', {
                                        'email': 'test@example.com',
                                        'firstName': 'test',
                                        'lastName': 'example'
                                        })
            response = views.mailinglist(request)
            response_json = json.loads(response.content)
            self.assertEqual(response.status_code, 500)
            self.assertFalse(response_json['success'])

    def test_password_reset_with_invalid_email(self):
        response = self.client.post('/accounts/password_reset/', {'email': 'nonexistant-email@example.com'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Email not found')

    def test_password_reset_valid_email(self):
        response = self.client.post('/accounts/password_reset/', {'email': self.test_user.email})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], True)

    def test_update_password_invalid_token(self):
        self.test_user.password_reset_token = '12345'
        some_time_in_the_future = timezone.now() + datetime.timedelta(1)
        self.test_user.password_token_expires = some_time_in_the_future
        self.test_user.save()

        invalid_token = '0987'
        response = self.client.post('/accounts/update_password/' + invalid_token + '/',
                                    {'password': 'new_password'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Invalid link')

    def test_update_password_expired_token(self):
        self.test_user.password_reset_token = '12345'
        some_time_in_the_past = timezone.now() - datetime.timedelta(1)
        self.test_user.password_token_expires = some_time_in_the_past
        self.test_user.save()

        response = self.client.post('/accounts/update_password/' + self.test_user.password_reset_token + '/',
                                    {'password': 'new_password'})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Invalid link')

    def test_update_password(self):
        self.test_user.password_reset_token = '12345'
        some_time_in_the_future = timezone.now() + datetime.timedelta(1)
        self.test_user.password_token_expires = some_time_in_the_future
        self.test_user.save()

        response = self.client.post('/accounts/update_password/' + self.test_user.password_reset_token + '/',
                                    {'password': 'new_password'})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], True)

    def create_test_user(self):
        test_user_fields = {'role': 0, 'role_other': 'role_other', 'city': 'city', 'country': 'country',
                            'email': 'email@example.com', 'firstName': 'firstName', 'has_image': False,
                            'hide_email': False, 'hide_number': True,
                            'institution': 'institution', 'lastName': 'lastName', 'password': 'password',
                            'phone_number': '123456789', 'state': 'state', 'title': 'title'}
        created_user = MyUser.objects.create_user(**test_user_fields)
        return created_user
