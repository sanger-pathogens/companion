#!/bin/bash
set -ev
echo "deb http://http.debian.net/debian testing main" > /etc/apt/sources.list.d/testing.list
echo "Package: *" >> /etc/apt/preferences
echo "Pin: release a=testing" >> /etc/apt/preferences
echo "Pin-Priority: -1" >> /etc/apt/preferences