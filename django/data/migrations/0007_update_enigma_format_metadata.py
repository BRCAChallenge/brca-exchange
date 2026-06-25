# -*- coding: utf-8 -*-


from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0006_update_release_notes'),
    ]

    operations = [
        # Reformat Enigma in sources
        migrations.RunSQL(
            """
                UPDATE data_release SET sources = replace(sources, 'Enigma', 'ENIGMA');
            """
        ),
    ]
