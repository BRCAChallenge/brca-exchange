# -*- coding: utf-8 -*-
# Generated by Django 1.9.4 on 2017-06-23 19:34


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0016_get_more_allele_frequency_data_from_exac'),
    ]

    operations = [
        migrations.AlterField(
            model_name='report',
            name='AA_Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_count_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_frequency_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Allele_number_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='EA_Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='report',
            name='Homozygous_count_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='AA_Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_count_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_frequency_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Allele_number_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='EA_Allele_Frequency_ESP',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_AFR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_AMR_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_EAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_FIN_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_NFE_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_OTH_ExAC',
            field=models.TextField(default=b'-'),
        ),
        migrations.AlterField(
            model_name='variant',
            name='Homozygous_count_SAS_ExAC',
            field=models.TextField(default=b'-'),
        ),
    ]
