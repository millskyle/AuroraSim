echo "Making HQ version..."
gs \
  -o "print_$1" \
  -sDEVICE=pdfwrite \
  -dDownsampleColorImages=true \
  -dDownsampleGrayImages=true \
  -dDownsampleMonoImages=true \
  -dColorImageResolution=2400 \
  -dGrayImageResolution=2400 \
  -dMonoImageResolution=2400 \
  -dColorImageDownsampleThreshold=1.0 \
  -dGrayImageDownsampleThreshold=1.0 \
  -dMonoImageDownsampleThreshold=1.0 \
   $1 > /dev/null
echo "Making screen version..."
gs \
  -o "screen_$1" \
  -sDEVICE=pdfwrite \
  -dDownsampleColorImages=true \
  -dDownsampleGrayImages=true \
  -dDownsampleMonoImages=true \
  -dColorImageResolution=300 \
  -dGrayImageResolution=300 \
  -dMonoImageResolution=300 \
  -dColorImageDownsampleThreshold=1.0 \
  -dGrayImageDownsampleThreshold=1.0 \
  -dMonoImageDownsampleThreshold=1.0 \
   "print_$1" > /dev/null

echo "Making web version..."
gs \
  -o "web_$1" \
  -sDEVICE=pdfwrite \
  -dDownsampleColorImages=true \
  -dDownsampleGrayImages=true \
  -dDownsampleMonoImages=true \
  -dColorImageResolution=120 \
  -dGrayImageResolution=120 \
  -dMonoImageResolution=120 \
  -dColorImageDownsampleThreshold=1.0 \
  -dGrayImageDownsampleThreshold=1.0 \
  -dMonoImageDownsampleThreshold=1.0 \
   "screen_$1" > /dev/null

















