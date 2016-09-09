import django
from django.contrib.auth.models import (
    BaseUserManager, AbstractBaseUser
)
from django.db import models
from django.utils import timezone

ROLE_OTHER = 0

class MyUserManager(BaseUserManager):
    def create_user(self, email, password, firstName="", lastName="", title="", role=ROLE_OTHER, role_other="", institution="",
                    city="", state="", country="", latitude="", longitude="", phone_number="", 
                    hide_number=False, hide_email=False, has_image=False,
                    is_admin=False, is_approved=False, admin_notifications=False):
        """
        Creates and saves a User with the given fields
        """
        if not email:
            raise ValueError('Users must have an email address')

        user = self.model(
            email=self.normalize_email(email),
            password=password,
            firstName=firstName,
            lastName=lastName,
            title=title,
            role=role,
            role_other=role_other,
            institution=institution,
            city=city,
            state=state,
            country=country,
            latitude=latitude,
            longitude=longitude,
            phone_number=phone_number,
            hide_number=hide_number,
            hide_email=hide_email,
            has_image=has_image,
            is_admin=is_admin,
            is_approved=is_approved,
            admin_notifications=admin_notifications,
        )

        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(is_admin=True, **kwargs)
        return user


class MyUser(AbstractBaseUser):
    ROLE_OTHER                  = ROLE_OTHER
    ROLE_COMMUNITY_MEMBER       = 1
    ROLE_CLINICAL_LAB_DIRECTOR  = 3
    ROLE_DIAGNOSTIC_LAB_STAFF   = 4
    ROLE_PRINCIPAL_INVESTIGATOR = 5
    ROLE_RESEARCHER             = 6
    ROLE_ADVOCACY_LEADER        = 7
    ROLE_ADVOCACY_MEMBER        = 8
    ROLE_GENETIC_COUNSELOR      = 9
    ROLE_CLINICAL_GENETICIST    = 10
    ROLE_CLINICIAN              = 11
    ROLE_DATA_PROVIDER          = 12

    email = models.EmailField(
        verbose_name='email address',
        max_length=255,
        unique=True,
    )
    firstName = models.TextField(blank=True)
    lastName = models.TextField(blank=True)
    title = models.TextField(blank=True)
    role = models.IntegerField(default=ROLE_OTHER)
    role_other = models.TextField(blank=True)
    institution = models.TextField(blank=True)
    city = models.TextField(blank=True)
    state = models.TextField(blank=True)
    country = models.TextField(blank=True)
    latitude = models.TextField(blank=True)
    longitude = models.TextField(blank=True)
    phone_number = models.TextField(max_length=30, blank=True)
    hide_number = models.BooleanField(default=False)
    hide_email = models.BooleanField(default=False)
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)
    is_approved = models.BooleanField(default=False)
    admin_notifications = models.BooleanField(default=False)
    has_image = models.BooleanField(default=False)

    activation_key = models.CharField(max_length=40, blank=True)

    password_reset_token = models.CharField(max_length=40, blank=True)
    password_token_expires = models.DateTimeField(default=django.utils.timezone.now)

    objects = MyUserManager()

    USERNAME_FIELD = 'email'

    def get_full_name(self):
        # The user is identified by their email address
        return self.email

    def get_short_name(self):
        # The user is identified by their email address
        return self.email

    def __str__(self):  # __unicode__ on Python 2
        return self.email

    def has_perm(self, perm, obj=None):
        "Does the user have a specific permission?"
        # Simplest possible answer: Yes, always
        return True

    def has_module_perms(self, app_label):
        "Does the user have permissions to view the app `app_label`?"
        # Simplest possible answer: Yes, always
        return True

    @property
    def is_staff(self):
        "Is the user a member of staff?"
        # Simplest possible answer: All admins are staff
        return self.is_admin

class MailingListEmail(models.Model):
    email = models.EmailField(
        verbose_name='email address',
        max_length=255,
        unique=True,
    )

