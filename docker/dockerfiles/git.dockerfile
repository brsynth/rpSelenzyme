FROM alpine

RUN apk -q update \
 && apk -q add git bash
