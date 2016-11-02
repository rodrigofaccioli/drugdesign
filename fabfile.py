#  Instalação do projeto:
#  
#  Realizar download do virtualenvwrapper
#  Criar virtualenv para o drugdesign
#  Instalar o requirements.txt
#  Rodar o arquivo fabfile.py

import os
import argparse
from fabric.api import run, sudo
from fabric.contrib.files import append, exists

#  path utils
CURRENT_PATH = os.path.dirname(os.path.abspath(__file__))
PROGRAMS_PATH = os.getenv("HOME") + '/programs/'

#  dependencies urls
spark_url = "http://d3kbcqa49mib13.cloudfront.net/\
            spark-1.6.2-bin-hadoop2.6.tgz"
mgltools_url = "http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/\
                mgltools_x86_64Linux2_1.5.6.tar.gz"
autodock_vina_url = "http://vina.scripps.edu/download/\
                    autodock_vina_1_1_2_linux_x86.tgz"


def get_params():
    parser = argparse.ArgumentParser(description="Script to deploy drugdesign")
    parser.add_argument("-h", action="store", dest="hosts", default=False,
                        required=False, help="e.g. 192.168.1.1, 192.168.0.1")
    parser.add_argument("-d", action="store", dest="dependencies",
                        default="dependencies.txt", required=False,
                        help="e.g. dependencies.txt model")
    parser.add_argument("-p", action="store", dest="packages",
                        default="packages.txt", required=False,
                        help="e.g. packages.txt model")
    parser.add_argument("--build_local", action="store", dest="build_local",
                        default=False, required=False, help="e.g. True")


#  Function to build configuration files
def build_config():
    return


#  Function to get OS
def get_os():
    return


#  Function to install OS packages
def install_packages():
    return


#  Function to download dependencies
def get_dependencies_by_url():
    return


#  Function to move files
def move_file():
    return


#  Function to unpack files tar -xf file.ext
def unpack_file():
    return


#  Function to install dependencies
def install_dependencies():
    return


#  Function to test dependencies
def test_dependencies():
    return


#  Function to compile C stuff
def install_drugdesign():
    return


#  Function to build drugdesign in local mode
def build_local():
    return


#  Function to build drugdesign in remote mode
def build_remote():
    return
