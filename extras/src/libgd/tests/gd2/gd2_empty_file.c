/* $Id: gd2_empty_file.c,v 1.1.1.1 2014/02/06 18:14:22 duncan Exp $ */
#include "gd.h"
#include <stdio.h>
#include <stdlib.h>
#include "gdtest.h"

int main()
{
	gdImagePtr im;
	FILE *fp;
	char path[1024];

	sprintf(path, "%s/gd2/empty.gd2", GDTEST_TOP_DIR);

	fp = fopen(path, "rb");
	if (!fp) {
		printf("failed, cannot open file (%s)\n", path);
		return 1;
	}

	im = gdImageCreateFromGd2(fp);
	fclose(fp);

	if (!im) {
		return 0;
	} else {
		gdImageDestroy(im);
		return 1;
	}
}
