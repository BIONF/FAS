#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

bashFile='.bashrc'
if [ "$sys" == "Darwin" ]; then
	shell=$(echo $SHELL)
	if [ $shell == "/bin/zsh" ]; then
    	bashFile='.zshrc'
	else
		bashFile='.bash_profile'
	fi
fi

### prepare folders
echo "-------------------------------------"
echo "Preparing folders..."
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR/..
CURRENT=$(pwd)

# create required folders
folders=(
  "CAST"
  "COILS2"
  "Pfam"
  "Pfam/Pfam-hmms"
  "Pfam/output_files"
  "SEG"
  "SignalP"
  "SMART"
  "TMHMM"
)

for i in "${folders[@]}"; do
  echo "$i"
  if [ ! -d $i ]; then mkdir $i; fi
done
echo "done!"

### download tools
echo "-------------------------------------"
echo "Downloading and installing annotation tools/databases:"
cd $CURRENT
if ! [ "$(ls -A $CURRENT/Pfam/Pfam-hmms)" ]; then
	if [ ! -f $CURRENT/data_HaMStR.tar ]; then
		echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
		wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
	else
		CHECKSUM=$(cksum data_HaMStR.tar)
		echo "Checksum: $CHECKSUM"
		if ! [ "$CHECKSUM" == "4100986910 5840435200 data_HaMStR.tar" ]; then
    		  rm $CURRENT/data_HaMStR.tar
    		  echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
      		  wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
    	fi
    fi

	if [ ! -f $CURRENT/data_HaMStR.tar ]; then
	  echo "File data_HaMStR.tar not found! Please try to download again from"
	  echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
	  exit
	fi

	CHECKSUM=$(cksum data_HaMStR.tar)
	if [ "$CHECKSUM" == "4100986910 5840435200 data_HaMStR.tar" ]; then
	  echo "Extracting archive data_HaMStR.tar"
	  tar xfv $CURRENT/data_HaMStR.tar
	  echo "Archive data_HaMStR.tar extracted into $CURRENT"
	  if [ ! -d $CURRENT/data_HaMStR ]; then
		echo "Directory $CURRENT/data_HaMStR not found!"
	  else
		printf "\nMoving Pfam ...\n---------------\n"
		rsync -rva data_HaMStR/Pfam/* $CURRENT/Pfam
		printf "\nMoving SMART ...\n----------------\n"
		rsync -rva data_HaMStR/SMART/* $CURRENT/SMART
		printf "\nMoving CAST ...\n---------------\n"
		rsync -rva data_HaMStR/CAST/* $CURRENT/CAST
		printf "\nMoving COILS ...\n----------------\n"
		rsync -rva data_HaMStR/COILS2/* $CURRENT/COILS2
		printf "\nMoving SEG ...\n--------------\n"
		rsync -rva data_HaMStR/SEG/* $CURRENT/SEG
		printf "\nMoving SignalP ...\n------------------\n"
		rsync -rva data_HaMStR/SignalP/* $CURRENT/SignalP
		printf "\nMoving TMHMM ...\n----------------\n"
		rsync -rva data_HaMStR/TMHMM/* $CURRENT/TMHMM
		rsync -rva data_HaMStR/README* $CURRENT/
		printf "\nRemoving duplicated data. Please wait.\n------------------------------------\n"
		rm -rf $CURRENT/data_HaMStR
  	    # rm $CURRENT/data_HaMStR.tar
		printf "\nDone! Data should be in place to run FAS.\n"
	  fi
	else
	  echo "Something went wrong with the download. Checksum does not match."
	  echo "Please try to download again from"
	  echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
	  echo "Please put it into $CURRENT folder and run this setup again!"
	  exit
	fi
fi

### add paths to bash profile file
echo "-------------------------------------"
echo "Adding path to ~/$bashFile"

if [ -z "$(grep PATH=$CURRENT ~/$bashFile)" ]; then
	echo "export PATH=$CURRENT:\$PATH" >> ~/$bashFile
fi

echo "Done! Now you can use FAS tool :-)"
