#! /bin/bash
mkdir -p "build_dir/tests"
mkdir -p "build_dir/examples"

cmake -S "tests" -B "build_dir/tests"
cmake --build "build_dir/tests"
cmake --build "build_dir/tests" -t test

cmake -S examples -B "build_dir/examples"
cmake --build "build_dir/examples"
