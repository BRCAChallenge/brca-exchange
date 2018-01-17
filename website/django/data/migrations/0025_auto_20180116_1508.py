# -*- coding: utf-8 -*-
# Generated and customized by zfisch on 2018-01-10 16:04

'''
Adds a name field which starts at 1 and is set in ascending
order based on the date field (when the release was generated).
'''

from __future__ import unicode_literals

from django.db import migrations, models
from data.utilities import set_release_name_defaults


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0024_auto_20171005_1906'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='datarelease',
            options={'ordering': ['date']},
        ),
        migrations.AddField(
            model_name='datarelease',
            name='name',
            field=models.PositiveIntegerField(null=True)
        ),
        migrations.RunPython(set_release_name_defaults),

        # make field non-nullable after defaults are set
        migrations.AlterField(
            model_name='datarelease',
            name='name',
            field=models.PositiveIntegerField(),
        ),
    ]
