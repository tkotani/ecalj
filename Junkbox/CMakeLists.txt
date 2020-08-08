CMakeLists.txt

cmake_minimum_required(VERSION 2.8)
# cmake 自体のバージョンを明示します。

enable_language(Fortran)
# cmake の名の通り、デフォルトでは C/C++ しか扱わないので、
# Fortran向けの設定を生成するように指定します。

project(simple-dependency Fortran)

set(EXECUTABLE abcd)
# CMakeの変数宣言です。実行体の名前を後から変更しやすいよう
# 最初から変数で定義したほうがいいと思います。

add_executable(${EXECUTABLE} # "abcd" のことです
  a.f90
  b.f90
  c.f90
  d.f90
  main.f90
  )
