# -*- coding: utf-8 -*-
# Generated by Django 1.11.10 on 2018-06-15 18:30
from django.db import migrations, models
import django_rdkit.models.fields


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('django_rdkit', '0001_setup'),
    ]

    operations = [
        migrations.CreateModel(
            name='BfpModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('bfp', django_rdkit.models.fields.BfpField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='CtabModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ctab', models.TextField(blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='MoleculeModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('molecule', django_rdkit.models.fields.MolField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='ReactionModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('rxn', django_rdkit.models.fields.RxnField()),
            ],
        ),
        migrations.CreateModel(
            name='SfpModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('sfp', django_rdkit.models.fields.SfpField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='SmartsModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('smarts', models.CharField(blank=True, max_length=2048)),
            ],
        ),
        migrations.CreateModel(
            name='SmilesModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('smiles', models.CharField(blank=True, max_length=2048)),
                ('molecule', django_rdkit.models.fields.MolField(null=True)),
            ],
        ),
    ]
