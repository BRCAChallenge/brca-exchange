ALLOWED_HOSTS = (
    'localhost:8080',
    'localhost'
)

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '&wy54&8++2d+%57ipbcuw#utl0f@!9sd6j*r1-cy@ihgz1%eae'
CAPTCHA_SECRET = '6LdwNBwTAAAAABi4FQqI_W8qzPiuDypc4re1DjuI'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

URL_FRONTEND = 'http://localhost:8080/'

EMAIL_HOST = ''
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''

# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'HOST': 'localhost',
        'NAME': 'storage.pg',
        'USER': 'postgres',
        'PASSWORD': 'postgres'
    }
}

