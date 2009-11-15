#! /bin/sh
./test_input_conversion_64_to_gmp 0 mat 17 || exit 1
./test_input_conversion_64_to_gmp 0 mat 1234567890123456 || exit 1
