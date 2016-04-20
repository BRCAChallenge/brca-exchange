# Create your tests here.
import datetime
import json

from django.test import TestCase
from django.utils import timezone

from .models import MyUser


class UserTestCase(TestCase):
    def test_password_reset_with_invalid_email(self):
        response = self.client.post('/accounts/password_reset/', {'email': 'nonexistant-email@example.com'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Email not found')

    def test_password_reset_valid_email(self):
        test_user = self.create_test_user()

        response = self.client.post('/accounts/password_reset/', {'email': test_user.email})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], True)

    def test_update_password_invalid_token(self):
        test_user = self.create_test_user()
        test_user.password_reset_token = '12345'
        some_time_in_the_future = timezone.now() + datetime.timedelta(1)
        test_user.password_token_expires = some_time_in_the_future
        test_user.save()

        invalid_token = '0987'
        response = self.client.post('/accounts/update_password/' + invalid_token + '/',
                                    {'password': 'new_password'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Invalid link')

    def test_update_password_expired_token(self):
        test_user = self.create_test_user()
        test_user.password_reset_token = '12345'
        some_time_in_the_past = timezone.now() - datetime.timedelta(1)
        test_user.password_token_expires = some_time_in_the_past
        test_user.save()

        response = self.client.post('/accounts/update_password/' + test_user.password_reset_token + '/',
                                    {'password': 'new_password'})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], False)
        self.assertEqual(json.loads(response.content)['error'], 'Invalid link')

    def test_update_password(self):
        test_user = self.create_test_user()
        test_user.password_reset_token = '12345'
        some_time_in_the_future = timezone.now() + datetime.timedelta(1)
        test_user.password_token_expires = some_time_in_the_future
        test_user.save()

        response = self.client.post('/accounts/update_password/' + test_user.password_reset_token + '/',
                                    {'password': 'new_password'})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['success'], True)

    def create_test_user(self):
        test_user_fields = {'affiliation': 'affiliation', 'city': 'city', 'comment': 'comment', 'country': 'country',
                            'email': 'email@example.com', 'firstName': 'firstName', 'has_image': False,
                            'hide_email': False, 'hide_number': True, 'include_me': True, 'email_me': True,
                            'institution': 'institution', 'lastName': 'lastName', 'password': 'password',
                            'phone_number': '123456789', 'state': 'state', 'title': 'title'}
        created_user = MyUser.objects.create_user(**test_user_fields)
        return created_user
