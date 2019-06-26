#!/bin/sh

sed -e 's/fdxrd_c/fdxr_42_c/g' -e 's/fdxox_c/fdxo_42_c/g' ./pre_update_icbi655_central.json > icbi655_central.json

