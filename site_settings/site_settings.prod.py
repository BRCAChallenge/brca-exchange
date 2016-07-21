
# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '' # enter secret key here
CAPTCHA_SECRET = '' # enter captcha secret here 

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False


URL_FRONTEND = 'http://brcaexchange.org/'

EMAIL_HOST = 'gmail.com'
EMAIL_HOST_USER = 'brcaexchange'
EMAIL_HOST_PASSWORD = '' # enter email password here

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'HOST': 'localhost',
        'NAME': 'production.pg',
        'USER': 'postgres',
        'PASSWORD': '' # enter database password here
    }
}

