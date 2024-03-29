# -*- coding: utf-8 -*-
# Generated by Django 1.9.13 on 2022-04-28 19:19
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0015_auto_20161011_1516'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mailinglistemail',
            name='email',
            field=models.EmailField(max_length=255, unique=True, verbose_name='email address'),
        ),
        migrations.AlterField(
            model_name='myuser',
            name='email',
            field=models.EmailField(max_length=255, unique=True, verbose_name='email address'),
        ),
    ]
