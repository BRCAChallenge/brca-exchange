from __future__ import unicode_literals

from django.db import migrations

def add_change_types(apps, schema_editor):
    ChangeType = apps.get_model("data", "ChangeType")
    ChangeType.objects.create(name="new")
    ChangeType.objects.create(name="deleted")
    ChangeType.objects.create(name="added_classification")
    ChangeType.objects.create(name="changed_classification")
    ChangeType.objects.create(name="added_information")
    ChangeType.objects.create(name="changed_information")

class Migration(migrations.Migration):

    dependencies = [
        ('data', '0009_change_type'),
    ]

    operations = [
        migrations.RunPython(add_change_types)
    ]
