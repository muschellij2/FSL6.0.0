#include <stdio.h>
#include <math.h>
#include "gd.h"
#include "gdtest.h"

#ifdef _MSC_VER
# define snprintf _snprintf
#endif

int gen_image(const char* filename, int idx, int reverse_x, int width, int height, int bgd)
{
	double gradient = height / (width*2.0);
	int offset = idx*width;
	int x1,y1,x2,y2, error = 0;

	gdImagePtr im = gdImageCreateTrueColor(width,height);
	if (bgd==1) {
		gdImageFilledRectangle(im,0,0,width-1,height-1, gdTrueColorAlpha(255, 255, 255, 0));
	} else {
		gdImageFilledRectangle(im,0,0,width-1,height-1, gdTrueColorAlpha(255, 255, 0, 0));
	}

	gdImageSetAntiAliased(im, gdTrueColorAlpha(0,0,0,0));

	/*
	test for potential segfault (introduced with AA improvements, fixed
	with the same patch - but I didn't notice it until later).*/
	gdImageLine(im,-1,-1,-1,-1,gdAntiAliased);

	x1 = floor(reverse_x * -width + 0.5);
	y1 = (offset-width) * gradient + 0.5;

	x2 = floor(reverse_x *  width*2 + 0.5);
	y2 = floor((offset+width*2) * gradient + 0.5);

	/* draw an AA line */
	gdImageLine(im, x1, y1, x2, y2, gdAntiAliased);

	gdImageLine(im, 0, im->sy - 1, im->sx, im->sy - 1, 0x40FF0000);

	if (!gdAssertImageEqualsToFile(filename, im)) {
		gdTestErrorMsg("gdAssertImageEqualsToFile failed: <%s>\n", filename);
		error = 1;
	}
	gdImageDestroy(im);
	return error;
}

int main()
{
	int error = 0;
	char path[1024];

	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_a_0_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,0,1,10,100, 1);

	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_a_1_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,1,1,10,100, 2);

	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_b_0_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,2,-1,10,100, 1);
	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_b_1_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,1,-1,10,100, 2);

	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_c_0_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,0,1,100,10, 1);
	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_c_1_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,1,1,100,10, 2);

	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_d_0_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,2,-1,100,10, 1);
	snprintf(path, sizeof(path)-1, "%s/gdimageline/gdimageline_aa_d_1_exp.png", GDTEST_TOP_DIR);
	error |= gen_image(path,1,-1,100,10, 2);

	return error;
}
