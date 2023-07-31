#!/bin/bash

coverage run -m pytest tests
coverage report -m
coverage html
