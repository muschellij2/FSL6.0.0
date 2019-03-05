/* $Id: bug00166.c,v 1.1.1.1 2014/02/06 18:14:22 duncan Exp $ */
#include "gd.h"
#include <stdio.h>
#include <stdlib.h>
#include "gdtest.h"

int
main(void)
{
	gdImagePtr im;
	char path[1024];
	int c, result;

	sprintf(path, "%s/xpm/bug00166.xpm", GDTEST_TOP_DIR);
	im = gdImageCreateFromXpm(path);
	if (!im) {
		return 2;
	}
	c = gdImageGetPixel(im, 1, 1);
	if (gdImageRed(im, c)      == 0xAA
	        && gdImageGreen(im, c) == 0xBB
	        && gdImageBlue(im, c)  == 0xCC) {
		result = 0;
	} else {
		result = 1;
	}
	gdImageDestroy(im);
	return result;
}
