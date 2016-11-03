#!/bin/bash/env python
# -*- coding: utf-8 -*-
import os
from fabric.api import *
from fabric.contrib.files import exists

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROGRAMS_PATH = os.getenv("HOME") + "/programs/"
SPARK_PATH = PROGRAMS_PATH + "spark/"
MGLTOOLS_PATH = PROGRAMS_PATH + "mgltools/"
AUTODOCKVINA_PATH = PROGRAMS_PATH + "autodock-vina/"

spark_url = "http://d3kbcqa49mib13.cloudfront.net/spark-1.6.2-bin-hadoop2.6.tgz"
mgltools_url = "http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz"
autodock_vina_url = "http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz"


def make_directory(directory):
    if not exists(directory):
        run("mkdir %s" % directory)


def get_dependency_by_url(url):
    run("wget %s" % url)


def move_file(file, path):
    run("mv %s %s" % (file, path))


def copy_file(file, path):
    run("cp -r %s %s" % (file, path))


def remove_file(file):
    run("rm -rf %s" % file)


def unpack_dependency(filename):
    run("tar -xf %s" % filename)


def install_packages():
    sudo("apt-get update")
    sudo("apt-get -y upgrade")
    sudo("apt-get -y install python-dev")
    sudo("apt-get -y install python-pip")
    sudo("apt-get -y install python-numpy")
    sudo("apt-get -y install python-matplotlib")
    sudo("apt-get -y install python-virtualenv")
    sudo("apt-get -y install pymol")
    sudo("pip install virtualenvwrapper")
    sudo("apt-get -y install zip")
    sudo("apt-get -y install tar")
    sudo("apt-get -y install openmpi-bin libopenmpi-dev")
    sudo("apt-get -y install cmake")
    sudo("apt-get -y install libffi-dev")
    sudo("apt-get -y install openjdk-8-jdk")


def install_dependencies():
    make_directory(PROGRAMS_PATH)

    with cd(PROGRAMS_PATH):
        get_dependency_by_url(spark_url)
        get_dependency_by_url(mgltools_url)
        get_dependency_by_url(autodock_vina_url)

        make_directory(SPARK_PATH)
        make_directory(MGLTOOLS_PATH)
        make_directory(AUTODOCKVINA_PATH)

        move_file("spark-1.6.2-bin-hadoop2.6.tgz", SPARK_PATH)
        move_file("mgltools_x86_64Linux2_1.5.6.tar.gz", MGLTOOLS_PATH)
        move_file("autodock_vina_1_1_2_linux_x86.tgz", AUTODOCKVINA_PATH)

        with cd(SPARK_PATH):
            unpack_dependency("spark-1.6.2-bin-hadoop2.6.tgz")

            with cd("spark-1.6.2-bin-hadoop2.6"):
                copy_file("*", SPARK_PATH)

            remove_file("spark-1.6.2-bin-hadoop2.6.tgz")
            remove_file("spark-1.6.2-bin-hadoop2.6")

        with cd(MGLTOOLS_PATH):
            unpack_dependency("mgltools_x86_64Linux2_1.5.6.tar.gz")

            with cd("mgltools_x86_64Linux2_1.5.6"):
                copy_file("*", MGLTOOLS_PATH)

            remove_file("mgltools_x86_64Linux2_1.5.6.tar.gz")
            remove_file("mgltools_x86_64Linux2_1.5.6")
            run("./install.sh -d %s" % MGLTOOLS_PATH)

        with cd(AUTODOCKVINA_PATH):
            unpack_dependency("autodock_vina_1_1_2_linux_x86.tgz")

            with cd("autodock_vina_1_1_2_linux_x86"):
                copy_file("*", AUTODOCKVINA_PATH)

            remove_file("autodock_vina_1_1_2_linux_x86.tgz")
            remove_file("autodock_vina_1_1_2_linux_x86")


def compile_drugdesign():
    with cd((CURRENT_DIR + "/virtualscreening/vina/mpi/")):
        make_directory("build")

        with cd("build"):
            run("cmake ..")
            run("make")

    with cd((CURRENT_DIR + "/virtualscreening/vina/detect_hbonds/")):
        make_directory("build")

        with cd("build"):
            run("cmake ..")
            run("make")


def build_drugdesign():
    install_packages()
    install_dependencies()
    compile_drugdesign()
