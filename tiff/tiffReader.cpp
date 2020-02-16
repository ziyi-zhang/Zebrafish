#include "tiffio.h"
int main(int argc, char* argv[])
{
    TIFF* tif = TIFFOpen(argv[1], "r");
    if (tif) {
	int dircount = 0;
	do {
	    dircount++;
	} while (TIFFReadDirectory(tif));
	printf("%d directories in %s\n", dircount, argv[1]);
	TIFFClose(tif);
    }
    return(0);
}
