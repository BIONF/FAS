# #!/bin/bash
# set -e

$PYTHON setup.py install --single-version-externally-managed --record=record.txt    # Python command to install the script.

#BINARY_HOME=$PREFIX/bin
#TESTPKG_HOME=$PREFIX/$PKG_NAME-$PKG_VERSION
#
## Copy source to the conda environment
#mkdir -p $TESTPKG_HOME
#cp -R $SRC_DIR/* $TESTPKG_HOME/
#
#echo "INSTALL FAS..."
#$PYTHON $TESTPKG_HOME/bin/setup.py install --single-version-externally-managed --record=record.txt
