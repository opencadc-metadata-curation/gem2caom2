#!/bin/bash

# IMAGE="bucket.canfar.net/gem2caom2"
IMAGE="gem_build_test"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
# docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run -m=7g --rm --name gem_run_builder -v ${PWD}:/usr/src/app/ ${IMAGE} gem_run || exit $?

date
exit 0