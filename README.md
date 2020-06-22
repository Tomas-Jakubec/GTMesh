[![pipeline status](https://mmg-gitlab.fjfi.cvut.cz/gitlab/jakubec/UnstructuredMesh/badges/master/pipeline.svg)](https://mmg-gitlab.fjfi.cvut.cz/gitlab/jakubec/UnstructuredMesh/commits/master)

# GTMesh

The GTMesh is a C++ library utilizing modern C++ paradigms as template metaprogramming
and type traits. The aim of GTMesh is to provide an implementation working with an unstructured mesh
of any dimension and topology. Furthermore, the library provides additional tools developed
during the development of UnstructuredMesh and its functionalities.

The tools developed as part of GTMesh:
- [unstructured mesh](src/GTMesh/UnstructuredMesh/) with simple and user friendly inferface
- user friendly, simple and generic [debugging tool](src/GTMesh/Debug/)
- [Traits](src/GTMesh/Traits/), a tool describing C++ data structures and
