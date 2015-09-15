#!/bin/sh

echo $$ > $1.pid
exec "$@"