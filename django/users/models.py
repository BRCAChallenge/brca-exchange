from django.contrib.auth.models import (
    BaseUserManager, AbstractBaseUser
)
from django.db import models


class MyUserManager(BaseUserManager):
    def create_user(self, email, password, firstName="", lastName="", title="", affiliation="", institution="",
                    city="", state="", country="",phone_number="", comment="", 
                    include_me=False, email_me=False, hide_number=False, hide_email=False, has_image=False,
                    is_admin=False, is_approved=False):
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
            affiliation=affiliation,
            institution=institution,
            city=city,
            state=state,
            country=country,
            phone_number=phone_number,
            comment=comment,
            include_me=include_me,
            email_me=email_me,
            hide_number=hide_number,
            hide_email=hide_email,
            has_image=has_image,
            is_admin=is_admin,
            is_approved=is_approved
        )

        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_superuser(self, **kwargs):
        user = self.create_user(is_admin=True, **kwargs)
        return user


class MyUser(AbstractBaseUser):
    email = models.EmailField(
        verbose_name='email address',
        max_length=255,
        unique=True,
    )
    firstName = models.TextField(blank=True)
    lastName = models.TextField(blank=True)
    title = models.TextField(blank=True)
    affiliation = models.TextField(blank=True)
    institution = models.TextField(blank=True)
    city = models.TextField(blank=True)
    state = models.TextField(blank=True)
    country = models.TextField(blank=True)
    phone_number = models.TextField(max_length=30, blank=True)
    hide_number = models.BooleanField(default=False)

    comment = models.TextField(blank=True)
    hide_email = models.BooleanField(default=False)

    include_me = models.BooleanField(default=True)
    email_me = models.BooleanField(default=True)
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)
    is_approved = models.BooleanField(default=False)
    has_image = models.BooleanField(default=False)

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
