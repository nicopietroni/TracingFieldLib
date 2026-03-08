#!/bin/bash

#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <subfolder_path>"
  exit 1
fi

SUBFOLDER="$1"
DONEFOLDER="$SUBFOLDER/done"

# Crea la sottocartella "done" se non esiste
mkdir -p "$DONEFOLDER"

# Cerca solo i .obj nella cartella principale (non ricorsivo)
find "$SUBFOLDER" -maxdepth 1 -type f -name "*.obj" | while read -r file; do
  echo "Processing: $file"
  ./quad_commandline "$file" -nofeatures -symmetry

  # Sposta il file nella sottocartella "done"
  BASENAME=$(basename "$file")
  mv "$file" "$DONEFOLDER/$BASENAME"
  echo "Moved $BASENAME to $DONEFOLDER"
done