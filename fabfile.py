#!/bin/bash/env python
# -*- coding: utf-8 -*-
import os
from fabric.api import *
from fabric.contrib.files import exists

PROGRAMS_PATH = os.getenv("HOME") + "/programs/"
SPARK_PATH = PROGRAMS_PATH + "spark/"
MGLTOOLS_PATH = PROGRAMS_PATH + "mgltools/"
AUTODOCKVINA_PATH = PROGRAMS_PATH + "autodock-vina/"
MPI_PATH = os.path.dirname(os.path.abspath(__file__)) + "/virtualscreening/vina/mpi/"
HBONDS_PATH = os.path.dirname(os.path.abspath(__file__)) + "/virtualscreening/vina/detect_hbonds/"

spark_url = "http://d3kbcqa49mib13.cloudfront.net/spark-1.6.2-bin-hadoop2.6.tgz"
mgltools_url = "http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz"
autodock_vina_url = "http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz"


def install_packages():
    local("sudo apt-get update")
    local("sudo apt-get -y upgrade")
    local("sudo apt-get -y install python-dev")
    local("sudo apt-get -y install python-pip")
    local("sudo apt-get -y install python-numpy")
    local("sudo apt-get -y install python-matplotlib")
    local("sudo apt-get -y install python-virtualenv")
    local("sudo apt-get -y install pymol")
    local("sudo pip install virtualenvwrapper")
    local("sudo apt-get -y install zip")
    local("sudo apt-get -y install tar")
    local("sudo apt-get -y install openmpi-bin libopenmpi-dev")
    local("sudo apt-get -y install cmake")
    local("sudo apt-get -y install libcurl14-gnutls-dev")
    local("sudo apt-get -y install libffi-dev")
    local("sudo apt-get -y install openjdk-8-jdk")
    local("sudo ln -s /usr/lib/libmpi_cxx.so.0.0.1 /usr/lib/libmpi_cxx.so.1")
    local("sudo ln -s /usr/lib/openmpi/lib/libmpi_cxx.so.0.0.1 /usr/lib/openmpi/lib/libmpi_cxx.so.1")
    local("sudo ln -s /usr/lib/libmpi.so.0 /usr/lib/libmpi.so.1")
    local("sudo ln -s /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libmpi.so.1")


def make_directory(directory):
    if not exists(directory):
        local("mkdir %s" % directory)


def get_dependency_by_url(url):
    local("wget %s" % url)


def move_file(file, path):
    local("mv %s %s" % (file, path))


def unpack_dependency(filename):
    local("tar -xf %s" % filename)


def install_dependencies():
    if not exists(PROGRAMS_PATH):
        make_directory(PROGRAMS_PATH)

    with lcd(PROGRAMS_PATH):
        get_dependency_by_url(spark_url)
        get_dependency_by_url(mgltools_url)
        get_dependency_by_url(autodock_vina_url)

        make_directory(SPARK_PATH)
        make_directory(MGLTOOLS_PATH)
        make_directory(AUTODOCKVINA_PATH)

        move_file("spark-1.6.2-bin-hadoop2.6.tgz", SPARK_PATH)
        move_file("mgltools_x86_64Linux2_1.5.6.tar.gz", MGLTOOLS_PATH)
        move_file("autodock_vina_1_1_2_linux_x86.tgz", AUTODOCKVINA_PATH)

        with lcd(SPARK_PATH):
            unpack_dependency("spark-1.6.2-bin-hadoop2.6.tgz")
            with cd("spark-1.6.2-bin-hadoop2.6"):
                move_file("*", "..")

        with lcd(MGLTOOLS_PATH):
            unpack_dependency("mgltools_x86_64Linux2_1.5.6.tar.gz")

            with cd("mgltools_x86_64Linux2_1.5.6"):
                move_file("*", "..")

            local("./install.sh -d %s" % MGLTOOLS_PATH)

        with lcd(AUTODOCKVINA_PATH):
            unpack_dependency("autodock_vina_1_1_2_linux_x86.tgz")

            with cd("autodock_vina_1_1_2_linux_x86"):
                move_file("*", "..")


def install_drugdesign():
    with lcd(MPI_PATH):
        make_directory("build")

        with lcd("build"):
            local("cmake ../")
            local("make")

    with lcd(HBONDS_PATH):
        make_directory("build")

        with lcd("build"):
            local("cmake ../")
            local("make")


def main():
    install_packages()
    install_dependencies()
    install_drugdesign()


if __name__ == "__main__":
    main()
