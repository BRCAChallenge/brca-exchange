from __future__ import unicode_literals

from django.db import migrations

def add_change_types(apps, schema_editor):
    ChangeType = apps.get_model("data", "ChangeType")
    ChangeType.objects.create(name="new")
    ChangeType.objects.create(name="deleted")
    ChangeType.objects.create(name="new_classification")
    ChangeType.objects.create(name="modified")

class Migration(migrations.Migration):

    dependencies = [
        ('data', '0009_change_type'),
    ]

    operations = [
        migrations.RunPython(add_change_types)
    ]
