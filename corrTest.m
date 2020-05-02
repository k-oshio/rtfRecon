//
//  quick test routines
//

#import <RecKit/RecKit.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

#import <NumKit/NumKit.h>

void        test0();
void		test1();	// corr
void		test2();	// shift bug chk
void		test3();	// fourier interpolation
void		test4();	// 180 pair in golden angle set
void		test5();	// ESPIRiT (coil profile)
void		test6();	// non-orthogonal ft (with mask)
void		test7();	// auto-seg
void		test8();	// SVD
void		test9();	// coil sensitivity (ft based)(2D)
void		test10();	// map scaling test
void		test11();	// coil sensitivity (ft based)(3D)
void		test12();	// optical flow
void		test13();	// 2d unwrap
void		test14();	// CSF (making single image block, Rinkan-2)
void		test15();	// Nav read
void		test16();	// CSF (PCA)
void		test17();	// CSF (making single image block, Rinkan-1000)
void		test18();	// 2D k-traj correction
void		test19();	// plot laplacian with shift scale
void		test20();	// test routine for maxIndex step
void		test21();	// combine with thres
void		test22();	// SVD-based coil profile est
void		test23();	// POCS-based coil profile est
void		test24();	// POCS-based coil profile est (using only k-space center)
void		test25();	// "Parallel imaging" for radial sampling
void		test26();	// shift est from central view (for Stk Star)
void		test27();	// create axial and sagittal from coronal scale images
						//	(to save time)

RecImage	*radRecon(RecImage *raw);

int
main()
{
    @autoreleasepool {
    //    test0();
    //	test1();
	//	test2();
	//	test3();
	//	test4();
	//	test5();	// 2D
	//	test6();
	//	test7();
	//	test8();
	//	test9();
	//	test10();
	//	test11();
	//	test12();
	//	test13();
	//	test14();
	//	test15();
	//	test16();
	//	test17();
	//	test18();
	//	test19();
	//	test20();
	//	test21();
	//	test22();
	//	test23();
	//	test24();
	//	test25();
		test26();
	//	test27();
    }
}

void
test0()
{
    int     i;
    float   x, y;

    for (i = 0; i < 100; i++) {
        x = (float)i / 100.0;
        y = sqrt(1 - x) + pow(x, 2.0/3.0);
        printf("%f %f\n", x, y);
    }
    for (i = 0; i < 100; i++) {
        x = -(float)i / 100.0;
        y = -sqrt(1 - x) + pow(-x, 2.0/3.0);
        printf("%f %f\n", x, y);
    }
    for (i = 0; i < 100; i++) {
        x = (float)i / 100.0;
        y = -sqrt(1 - x) + pow(x, 2.0/3.0);
        printf("%f %f\n", x, y);
    }
    for (i = 0; i < 100; i++) {
        x = -(float)i / 100.0;
        y = -sqrt(1 - x) + pow(-x, 2.0/3.0);
        printf("%f %f\n", x, y);
    }
}

void
test1()	// corr
{
	RecImage    *pw, *pwr;
	RecImage    *roi, *roir;
	RecImage	*corr;
	NSPoint		pt;
	int			xdim = 32;
	int			ydim = 32;
	int			xpos = 112;	// [0..256]
	int			ypos = 48;	// [0..64]

	pw  = [RecImage imageWithKOImage:@"pw.img"];
	pwr = [RecImage imageWithKOImage:@"pwr2.img"];

	if (0) {
		[pw gauss2DHP:0.1];
		[pwr gauss2DHP:0.1];
		[pw saveAsKOImage:@"pw-hp.img"];
		[pwr saveAsKOImage:@"pwr2-hp.img"];
		exit(0);
	}

	roi  = [RecImage imageOfType:RECIMAGE_REAL xDim:xdim yDim:ydim];
	roir = [RecImage imageOfType:RECIMAGE_REAL xDim:xdim yDim:ydim];
	corr = [RecImage imageOfType:RECIMAGE_REAL xDim:xdim yDim:ydim];
	[roir copyLoopsOf:roi];
	[corr copyLoopsOf:roi];

	takeROI2d(pw,  roi,  xpos, ypos, 2);
	takeROI2d(pwr, roir, xpos, ypos, 2);
	[roi saveAsKOImage:@"roi1.img"];
	corr = [roi xyCorrelationWith:roir width:0.2 triFilt:YES];
	[corr saveAsKOImage:@"corr_64.img"];
	pt = Rec_find_peak2([corr data], xdim, ydim);
	printf("sft = (%f / %f)\n", pt.x, pt.y);
}

RecImage *
radRecon(RecImage *raw)
{
	RecImage	*traj, *img;
	RecGridder	*grid;
	RecLoop		*ky, *kx;

	traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:[raw yLoop], [raw xLoop], nil];
	kx = [RecLoop loopWithDataLength:256];
	ky = [RecLoop loopWithDataLength:256];
	img  = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[raw zLoop], ky, kx, nil];
	[traj initRadialTraj];          // for gridding
	grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
	[grid grid2d:raw to:img];

	return img;
}

void
dump_buf(float *p, float *q, int len)
{
	int		i;
	for (i = 0; i < len; i++) {
		printf("%d %f %f\n", i, p[i], q[i]);
	}
}

float
sinc(float th)
{
	if (th == 0) {
		return 1.0;
	} else {
		return sin(th) / th;
	}
}

// wiki: trigonometric interpolation
void
triginterp(float *yout, int n1, float *y, int n2)	// interval: 0::2PI
{
	int		i, k;
	float	dth1 = 2 * M_PI / n1;
	
	for (i = 0; i < n1; i++) {
		yout[i] = 0;
	}
int N = 4;
	for (k = 0; k < N; k++) {
		for (i = 0; i < n1; i++) {
			yout[i] += cos(dth1 * i * k) / N;
		}
	}
}

void
test4()	// 180 pair in golden angle set
{
	RecImage	*gaTab;
	RecLoop		*pe;
	int			i, j, nView;
	float		*tab;
	float		a1, a2, mn, dif;
	int			ix;

	pe = [RecLoop loopWithDataLength:256];
	gaTab = Rec_goldenAngleTab(pe, 8);
	nView = [gaTab xDim];
	tab = [gaTab data];

//	sort_mg(tab, nView);

//	for (i = 0; i < nView; i++) {
//		printf("%d %f\n", i, tab[i]);
//	}

	for (i = 0; i < nView; i++) {
		a1 = tab[i];
		mn = 1.0;
		ix = 0;
		for (j = 0; j < nView; j++) {
			if (j == i) continue;
			a2 = tab[j];
			dif = fabs(fabs(a1 - a2) - M_PI);
			if (dif < mn) {
				mn = dif;
				ix = j;
			}
		}
		printf("%d %d\n", i, ix);
	}
}

// ESPIRiT
void
test5()	// coil map est without ref
{
	RecImage		*img, *tmp_img, *raw, *map, *eval, *ref;
	RecImage		*g_im, *v_im;
	RecLoop			*ch, *sl;
	int				xDim, yDim, nCh;
	float			*p1, *q1;
	float			*p2, *q2;
	float			*p3, *q3;
	float			*p4, *q4;
	float			*p5;
	float			mx;
	int				x, y, c, i, j, ii, jj, ix1, ix2;
	RecImage		*tmp_im;
	Num_svd_result	*sres;
	Num_mat			*A;
	int				D = 64;	// 24	// size of calibration area
	int				d = 8;	// 6	// size of calibration kernel
	int				Dd = D - d + 1;	// nrow of calibration matrix A (nc = d x d x nCch)
	int				nEv = 100;		// num of singular values retained
	float			evThres = 1.0e-2;
	int				ftDim = 64;
	int				nr, nc;

// convert KO_IMAGE
if (0) {
	img     = [RecImage imageWithKOImage:@"/Users/oshio/Math/SMASH/spine_coil_phantom/IMG.re"];
	tmp_img = [RecImage imageWithKOImage:@"/Users/oshio/Math/SMASH/spine_coil_phantom/IMG.im"];
	[img makeComplexWithIm:tmp_img];
	[img saveAsKOImage:@"/Users/oshio/Math/SMASH/spine_coil_phantom/IMG.cpx"];

exit(0);
}

// rtf3d volume
if (1) {
	system("rm *.img");
//	system("rm IMG_*");
	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [img zLoop];
	img = [img sliceAtIndex:10 forLoop:sl];
	[img saveAsKOImage:@"img.img"];

	xDim = [img xDim];
	yDim = [img yDim];
	nCh = [ch dataLength];
}

// spine (single slice phantom image
if (0) {
	img = [RecImage imageWithKOImage:@"/Users/oshio/Math/SMASH/spine_coil_phantom/IMG.cpx"];
	xDim = [img xDim];
	yDim = [img yDim];
	ch = [img zLoop];
	nCh = [ch dataLength];
[img fft2d:REC_INVERSE];
[img takeEvenLines];
[img fft2d:REC_FORWARD];
}

// === single slice ====
	raw = [img copy];
	[raw fft2d:REC_INVERSE];
	tmp_img = [img combineForLoop:ch];
	[img saveAsKOImage:@"img.img"];
	[tmp_img saveAsKOImage:@"img_c.img"];
//	ftDim = [img xDim];

//	raw = [raw sliceAtIndex:32 forLoop:sl];
	[raw saveAsKOImage:@"raw.img"];
	[raw crop:[raw xLoop] to:D];
	[raw crop:[raw yLoop] to:D];
	[raw saveAsKOImage:@"IMG_cal"];

// fill A matrix
	nc = d * d * nCh;
	nr = Dd * Dd;
	tmp_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nc yDim:nr];
	p1 = [raw real];	// src
	q1 = [raw imag];
	p2 = [tmp_im real];	// dst
	q2 = [tmp_im imag];

// calibration matrix A (tmp_im)
	for (y = 0; y < Dd; y++) {
		for (x = 0; x < Dd; x++) {
			for (c = 0; c < nCh; c++) {
				for (ii = 0; ii < d; ii++) {
					for (jj = 0; jj < d; jj++) {
						ix1 = c * D*D + (y + ii) * D + x + jj;				// src
						ix2 = (y * Dd + x) * nc + c * d*d + ii * d + jj;	// dst
						p2[ix2] = p1[ix1];
						q2[ix2] = q1[ix1];
					}
				}
			}
		}
	}
	[tmp_im saveAsKOImage:@"IMG_A.img"];
	A = Num_im_to_m(tmp_im);

// SVD of A
	sres = Num_svd(A);
	Num_free_mat(A);
	saveAsKOImage(sres->Vt, @"IMG_Vt.img");
	saveAsKOImage(sres->U, @"IMG_U.img");

// singular value
	for (i = 0; i < sres->s->n; i++) {
		printf("%d %f\n", i, sres->s->data[i]);
	}
	mx = 0;
	for (i = 0; i < sres->s->n; i++) {
		mx = Rec_max(mx, sres->s->data[i]);
	}
	mx *= evThres;
	nEv = sres->s->n;
	for (i = 0; i < sres->s->n; i++) {
		if (sres->s->data[i] < mx) {
			nEv = i;
			break;
		}
	}
	printf("nEv = %d\n", nEv);
	
// Vt -> kernel image
	tmp_im = Num_m_to_im(sres->Vt);
	Num_free_svd_result(sres);

	v_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:d yDim:d zDim:nCh * nEv];
	p1 = [tmp_im real];	// src
	q1 = [tmp_im imag];
	p2 = [v_im real];	// dst
	q2 = [v_im imag];

	for (i = 0; i < d * d * nCh * nEv; i++) {
		p2[i] = p1[i];
		q2[i] = q1[i];
	}
//	[v_im saveAsKOImage:@"IMG_kern_f.img"];
	[v_im zeroFill:[v_im xLoop] to:ftDim];
	[v_im zeroFill:[v_im yLoop] to:ftDim];
	[v_im fft2d:REC_FORWARD];
	[v_im saveAsKOImage:@"IMG_kern.img"];

	g_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nCh yDim:nEv];
	tmp_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nCh yDim:nCh];
	map = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:ftDim yDim:ftDim zDim:nCh * nCh];
	eval = [RecImage imageOfType:RECIMAGE_REAL xDim:ftDim yDim:ftDim zDim:nCh];
	A = Num_im_to_m(g_im);
	sres = Num_new_svd_result(A);
	p1 = [v_im real];	// Vt of A
	q1 = [v_im imag];
	p2 = [g_im real];	// G
	q2 = [g_im imag];
	p3 = [tmp_im real];	// Vt of G
	q3 = [tmp_im imag];
	p4 = [map real];
	q4 = [map imag];
	p5 = [eval real];

	for (i = 0; i < ftDim; i++) {
		for (j = 0; j < ftDim; j++) {
			// make g_im
			for (ii = 0; ii < nEv; ii++) {
				for (jj = 0; jj < nCh; jj++) {
					// v_im : x:j, y:i, z:ii*nCh + jj
					// g_im : x:jj, y:ii
					p2[ii * nCh + jj] = p1[((ii * nCh + jj) * ftDim + i) * ftDim + j];
					q2[ii * nCh + jj] = q1[((ii * nCh + jj) * ftDim + i) * ftDim + j];
				}
			}
			// svd
			Num_copy_im_to_m(A, g_im);
			Num_svd_ref(A, sres->U, sres->s, sres->Vt);
			Num_copy_m_to_im(tmp_im, sres->Vt);
			// copy to map / eval
			for (ii = 0; ii < nCh; ii++) {		// eVal
				p5[(ii * ftDim + i) * ftDim + j] = sres->s->data[ii];	// OK
				for (jj = 0; jj < nCh; jj++) {	// ch
					p4[((ii * nCh + jj) * ftDim + i) * ftDim + j] = p3[ii * nCh + jj];
					q4[((ii * nCh + jj) * ftDim + i) * ftDim + j] = q3[ii * nCh + jj];
				}
			}
		}
	}
	// phase correction
	for (ii = 0; ii < nCh; ii++) {
		ref = [map sliceAtIndex:ii * nCh];
		[ref toUnitImage];
		[ref conjugate];
		for (jj = 0; jj < nCh; jj++) {
			tmp_im = [map sliceAtIndex:ii * nCh + jj];
			[tmp_im multByImage:ref];
			[map copySlice:tmp_im atIndex:ii * nCh + jj];
		}
	}
/*
	// interp to original resolution (not equivalent ???)
	[map fft2d:REC_INVERSE];
	[map zeroFill:[map xLoop] to:256];
	[map zeroFill:[map yLoop] to:256];
	[map fft2d:REC_FORWARD];
*/

	[map saveAsKOImage:@"IMG_map"];
	[eval saveAsKOImage:@"IMG_eval"];

	Num_free_mat(A);
	Num_free_svd_result(sres);
}

void
dump_p(float *p, int n)
{
	int	i;

	for (i = 0; i < n; i++) {
		printf("%d %f\n", i, p[i]);
	}
}

void
test6()	// non-orthogonal ft (with mask)
{
	RecImage	*img, *mask, *est, *err;
	int			i, n = 256;
	float		th, *p;
	float		w = 0.02;

	img = [RecImage imageOfType:RECIMAGE_REAL xDim:n];
	mask = [RecImage imageWithImage:img];

// test mask
	p = [mask data];
	for (i = 0; i < n; i++) {
		if (i > 50 && i < 190) {
			p[i] = 1.0;
		}
	}
//	dump_p(p, n);

// test data
	p = [img data];
	for (i = 0; i < n; i++) {
		th = (float)i * 2 * M_PI / n;
		p[i] = sin(th) + cos(th * 2) + sin(th * 20) * 0.2;
	}
	[img maskWithImage:mask];
//	dump_p(p, n);

// initial est
	est = [img copy];
	[est gauss1DLP:w forLoop:[est xLoop]];
//	dump_p([est data], n);

	for (i = 0; i < 2; i++) {
		err = [img copy];
		[err subImage:est];
		[err maskWithImage:mask];
		[err gauss1DLP:w forLoop:[err xLoop]];
	//	dump_p([err data], n);

		[est addImage:err];
	}
	[est maskWithImage:mask];
	dump_p([est data], n);
}

void
test7()
{
	RecImage	*img, *radMap, *rad;
	RecLoop		*ch, *sl;
	float		*p, v, w, sx, sy;
	int			i, j, n;

printf("start\n");
	system("rm *.img");
//	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav2.recimg" relativePath:YES];

	ch = [RecLoop findLoop:@"Channel"];

	img = [img combineForLoop:ch];
	img = [img avgForLoop:[img zLoop]];
	[img saveAsKOImage:@"img.img"];
	n = [img xDim];	// == yDim

	// find c.o.g.
	p = [img data];
	sx = sy = w = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			v = p[i * n + j];
			sx += j * v;
			sy += i * v;
			w += v;
		}
	}
	sx /= w;
	sy /= w;
	printf("%f %f\n", sx, sy);

	// to radial
//	rad = [img copy];
//	radMap = [img mapForPolarWithRMin:1 logR:NO];
//	[rad resample:img withMap:radMap];

	rad = [img toPolarWithNTheta:256 nRad:128 rMin:1 logR:NO];
	[rad saveAsKOImage:@"img_polar.img"];
}

void
test8()
{
	RecImage	*img, *raw;
	RecImage	*traj, *radMap;
	RecImage	*gaTab;
	RecImage	*gaSort;
	RecLoop		*ch, *sl, *pe, *rd, *kx, *ky;
	int			recDim		= 256;
	int			nRef		= 8;
	int			echoStart	= 112;

	Num_mat			*A, *U, *Vt;
	Num_vec			*s;
	int				nr, nc, n;

	printf("start\n");
	system("rm *.img");
	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-12-10/Run54033.5212.10-GoldenAngle"];	// g.a.

	rd = [RecLoop findLoop:@"Read"];
	pe = [RecLoop findLoop:@"Phase"];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [RecLoop findLoop:@"Slice"];
	kx = [RecLoop loopWithName:@"kx" dataLength:recDim];
	ky = [RecLoop loopWithName:@"ky" dataLength:recDim];
	
	sl = [raw cropToPo2:sl];		// implement dft (for small len) later ###

	gaTab = Rec_goldenAngleTab(pe, nRef);
	gaSort = Rec_goldenAngleSort(gaTab);
	
	[raw reorder:pe withTab:gaSort];
	img = [raw avgForLoop:ch];
	[img fft1d:sl direction:REC_FORWARD];
	[img saveAsKOImage:@"IMG_raw"];

	[img fft1d:rd direction:REC_FORWARD];
	[img saveAsKOImage:@"IMG_sin"];

//if (0) {
//RecImage *uimg, *simg, *vtimg;
//[img takeRealPart];
//[img recSVD_U:&uimg S:&simg Vt:&vtimg];
//[uimg saveAsKOImage:@"IMG_uimg"];
//[simg saveAsKOImage:@"IMG_simg"];
//[vtimg saveAsKOImage:@"IMG_vtimg"];
//exit(0);
//}

	img = [img sliceAtIndex:32];
	nr = [img yDim];
	nc = [img xDim];
	n = MIN(nr, nc);
	A = Num_im_to_m(img);
	// make U, s, Vt
	U = Num_new_cmat(nr, n);
	Vt = Num_new_cmat(n, nc);
	s = Num_new_cvec(n);
	// call LAPACK

//	Num_svd(A, U, s, Vt);	// 5.5 sec


	img = Num_v_to_im(s);
	[img saveAsKOImage:@"IMG_s"];

	img = Num_m_to_im(U);
	[img trans];
	[img saveAsKOImage:@"IMG_U"];

	img = Num_m_to_im(Vt);
	[img saveAsKOImage:@"IMG_Vt"];


//	[traj initRadialTrajWithSortTab:gaSort startPoint:echoStart];
//	[radMap  initRadialMapWithTab:gaSort];

}

void
test9()
{
	RecImage	*img, *tmp_img, *raw, *map, *sos;
	int			xDim, yDim, nCh;
	RecLoop		*ch, *sl;

printf("start\n");
	system("rm *.img");
	system("rm IMG_*.*");
	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [img zLoop];
	img = [img sliceAtIndex:32 forLoop:sl];
	[img saveAsKOImage:@"img.img"];

	xDim = [img xDim];
	yDim = [img yDim];
	nCh = [ch dataLength];

	raw = [img copy];

	// make test image
	if (0) {
		float		*p, *q;
		float		r, th, x, y;
		int			i, j, k, ix;
		RecNoise	*nz = [RecNoise noise];
		float		nlevel = 0.1; //0.1;	// 0.2 : ok, 0.3: ng

		p = [raw real];
		q = [raw imag];
		for (k = 0; k < nCh; k++) {
			for (i = 0; i < yDim; i++) {
				y = i - yDim/2;
				for (j = 0; j < xDim; j++) {
					x = j - xDim/2;
					r = sqrt(x*x + y*y);
					th = 0.005 * r * r;
					ix = k * xDim * yDim + i * xDim + j;
					if (r > 100) {
						p[ix] = [nz nrml] * nlevel;
						q[ix] = [nz nrml] * nlevel;
					} else {
						p[ix] = cos(th) + [nz nrml] * nlevel;
						q[ix] = sin(th) + [nz nrml] * nlevel;
					}
				}
			}
		}
		[raw saveAsKOImage:@"IMG_test.img"];
	}

//	[raw phaseUnwrap];

	[raw saveAsKOImage:@"IMG_lpf"];
	sos = [raw combineForLoop:[raw zLoop]];
	[sos saveAsKOImage:@"IMG_sos"];
	[raw divImage:sos];
	[raw saveAsKOImage:@"IMG_norm"];
}

RecImage *
mksub(RecImage *img, int dx, int dy, int xc, int yc, int slc)
{
    Rec2DRef    *ref;
	int			x0 = xc - dx/2;
	int			y0 = yc - dy/2;

    ref = [Rec2DRef refForImage:img];
    [ref setX:x0 y:y0 nX:dx nY:dy];
	[ref setZ:slc];
    return [ref makeImage];
}

//	offs : ok for odd sz, wrong for even sz
//	scale: wrong
//
void
test10()
{
	RecImage	*low, *hi;
	int			i, j;
	int			sz;
	float		*p;

	system("rm ../test_img/test10*.*");
//
	sz = 16;
	low = [RecImage imageOfType:RECIMAGE_REAL xDim:sz yDim:sz];
	p = [low data];

	for (i = 0; i < sz; i++) {
		if (i % 2 == 1) continue;
		for (j = 0; j < sz; j++) {
			p[i * sz + j] = j;
		}
	}
	[low saveAsKOImage:@"../test_img/test10-1.img"];
	hi = [low copy];
	hi = [hi scale1dLoop:[hi xLoop] by:256.0/sz to:256];	// sft is wrong ###
	hi = [hi scale1dLoop:[hi yLoop] by:256.0/sz to:256];
	[hi saveAsKOImage:@"../test_img/test10-2.img"];
	[hi gauss2DHP:0.2];
	[hi saveAsKOImage:@"../test_img/test10-3.img"];

	sz = 17;
	low = [RecImage imageOfType:RECIMAGE_REAL xDim:sz yDim:sz];
	p = [low data];

	for (i = 0; i < sz; i++) {
		if (i % 2 == 1) continue;
		for (j = 0; j < sz; j++) {
			p[i * sz + j] = j;
		}
	}
	[low saveAsKOImage:@"../test_img/test10-4.img"];
	hi = [low copy];
	hi = [hi scale1dLoop:[hi xLoop] by:256.0/sz to:256];	// sft is wrong ###
	hi = [hi scale1dLoop:[hi yLoop] by:256.0/sz to:256];
	[hi saveAsKOImage:@"../test_img/test10-5.img"];
}

void
test11()
{
	RecImage	*img, *mg, *raw, *sos, *lowres, *mask;
	int			xDim, yDim, zDim, nCh;
	int			res = 16;
	RecLoop		*ch, *sl;

printf("start\n");
	system("rm *.img");
	system("rm IMG_*.*");
	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [img zLoop];
	[img saveAsKOImage:@"img.img"];

	xDim = [img xDim];
	yDim = [img yDim];
	zDim = [img zDim];
	nCh = [ch dataLength];

	raw = [img copy];

	lowres = [raw copy];
	[lowres fermiWithRx:0.5 ry:0.5 d:0.1 x:0 y:0 invert:NO half:NO];
	[lowres fft3d:REC_INVERSE];
	[lowres crop:[lowres xLoop] to:res];
	[lowres crop:[lowres yLoop] to:res];
	[lowres crop:[lowres zLoop] to:res];
	[lowres saveAsKOImage:@"IMG_lowres_k.img"];
	[lowres fft3d:REC_FORWARD];
	[lowres saveAsKOImage:@"IMG_lowres.img"];
	mask = [lowres copy];
	[mask magnitude];
	[mask thresAt:0.05];
	[mask saveAsKOImage:@"IMG_mask.img"];
//mask = [mask sumForLoop:ch];
//[mask saveAsKOImage:@"IMG_mask_total.img"];
	

exit(0);


	[raw phaseUnwrap3D];	// imput is complex image
[raw saveAsKOImage:@"IMG_uwp.img"];
exit(0);

	mg = [img copy];
	[mg gauss3DLP:0.3];
	[mg magnitude];
	[mg gauss3DLP:0.03];

//	[raw phaseUnwrap3D];

	[raw saveAsKOImage:@"IMG_lpf"];
	[raw makeComplex];
	sos = [raw combineForLoop:ch];
	[sos saveAsKOImage:@"IMG_sos"];
	[raw divImage:sos];
	[raw saveAsKOImage:@"IMG_norm"];
}

void
test12()
{
	RecImage	*img1, *img2;
	RecImage	*sft;

	img1 = [RecImage imageWithKOImage:@"optical/imgb0_coro.img"];
	img2 = [RecImage imageWithKOImage:@"optical/imgb1_coro.img"];

	sft = [img1 opticalFlow3dWithRef:img2];
}

// CSF (phase)
void
test13()
{
	RecImage	*img;

	img = [RecImage imageWithKOImage:@"optical/b500-test-cpx-0"];

//	[img phaseUnwrap];
	[img saveAsKOImage:@"optical/b500-test-cpx-mg-phs"];
}

// CSF (mag)
void
test14()
{
	RecImage			*img, *img_m;
	NSArray				*imgLps;
	NSMutableArray		*loops;
	RecLoop				*cLoop;	// cardiac phase
	RecLoop				*avLoop;
	NSString			*base = @"../toshiba_images/DWI-rinkan-2";
	int					bval = 1000;
	NSString			*path;
	int					i;

// read images
	i = 0;
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-%d-cpx", base, bval, bval, i];
	img = [RecImage imageWithKOImage:path];
	avLoop = [img zLoop];
	cLoop = [RecLoop loopWithDataLength:7];
	imgLps = [img loops];
	loops = [NSMutableArray arrayWithArray:imgLps];
	[loops insertObject:cLoop atIndex:0];
	img_m = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loops];
	[img_m copySlice:img atIndex:0 forLoop:cLoop];

	for (i = 1; i < 7; i++) {
		printf("%d\n", i);
		path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-%d-cpx", base, bval, bval, i];
		img = [RecImage imageWithKOImage:path];
		[img setLoops:imgLps];
		[img_m copySlice:img atIndex:i forLoop:cLoop];
	}
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d.img", base, bval, bval];
	[img_m saveAsKOImage:path];
}

void
test15()
{
	FILE		*fp;
	int			i, isft;
	char		*buf, *sts;
	float		*sft;
	float		*hst;
	float		mn, mx;
	int			nbin = 40;
	int			nview;

	buf = (char *)malloc(256);
//	fp = fopen("../toshiba_images/2017-10-11/edge.txt", "r");
	fp = fopen("../toshiba_images/2017-10-24/edge10.txt", "r");
	if (fp == NULL) {
		printf("edge file not found\n");
		exit(0);
	}

	for (i = 0; ; i++) {
		sts = fgets(buf, 256, fp);
		if (sts == NULL) break;
	}
	nview = i;
	sft = (float *)malloc(sizeof(float) * nview);
	rewind(fp);
	for (i = 0; i < nview; i++) {
        sts = fgets(buf, 256, fp);
		if (sts == NULL) break;
		sscanf(buf, "[new Thread] - correctEdgePoint..%d", &isft);
		sft[i] = isft * 0.46875;
		printf("%d %f\n", i, sft[i]);
	}
	fclose(fp);

	mx = mn = sft[0];
	for (i = 0; i < nview; i++) {
		if (mn > sft[i]) mn = sft[i];
		if (mx < sft[i]) mx = sft[i];
	}
	printf("mn = %f, mx = %f\n", mn, mx);

	hst = (float *)malloc(sizeof(float) * nbin);
	histogram(hst, nbin, sft, nview, mn, mx * 1.0);

	for (i = 0; i < nbin; i++) {
		printf("%d %f\n", i, hst[i]);
	}
	free(hst);
	free(buf);
}

void
test16()
{
//	NSString			*base = @"../toshiba_images/DWI-rinkan-2";
	NSString			*base = @"../toshiba_images/DWI-rinkan-1024";
	int					bval = 100;
	int					slc = 1;
	NSString			*path;
	RecImage			*img, *mask;
	RecLoop				*tLoop;

	Num_svd_result		*sres;
	Num_ica_result		*ires;
	Num_mat				*A;
	RecImage			*a_im, *u_im, *vt_im;
	RecImage			*avg;
	int					np, n;
	BOOL				mag = NO;

	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-%d.img", base, bval, bval, slc];
//	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-%d-p.img", base, bval, bval, slc];
	img = [RecImage imageWithKOImage:path];
	if (img == nil) exit(0);

	if (mag) {
		[img magnitude];
	} else {
		// make mask
		mask = [img copy];
		[mask magnitude];
		[mask thresAt:0.1];
//[mask saveAsKOImage:@"IMG_unwrap.img"];
//exit(0);
		[img phase];
	//	[img multByImage:mask];
	//	[img unwrap2d];
	}
[img saveAsKOImage:@"IMG_unwrap.img"];
//exit(0);

	tLoop = [img zLoop];
	avg = [img avgForLoop:tLoop];
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_avg.img", base, bval, bval];
	[avg saveAsKOImage:path];
	// PCA
	[img subImage:avg];
	np = [img xDim] * [img yDim];
	n = [tLoop dataLength];
	a_im = [RecImage imageOfType:RECIMAGE_REAL xDim:np yDim:n];
//	a_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:np yDim:n];
	[a_im copyImageData:img];

	A = Num_im_to_m(a_im);
	sres = Num_svd(A);
	u_im = Num_m_to_im(sres->U);
	[u_im trans];
//	Num_scale_rows(sres->Vt, sres->s);
	vt_im = Num_m_to_im(sres->Vt);
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_U.img", base, bval, bval];
	[u_im saveAsKOImage:path];
	[img copyImageData:vt_im];
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_E.img", base, bval, bval];
	[img saveAsKOImage:path];

	int ncomp = 10;
	ires = Num_ica(A, ncomp);
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_out.img", base, bval, bval];
	saveAsKOImage(ires->WX, path);
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_W.img", base, bval, bval];
	saveAsKOImage(ires->W, path);
	vt_im = Num_m_to_im(ires->WK);
	img = [RecImage imageOfType:RECIMAGE_REAL xDim:[avg xDim] yDim:[avg yDim] zDim:ncomp];
	[img copyImageData:vt_im];
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d_map.img", base, bval, bval];
	[img saveAsKOImage:path];

	Num_free_svd_result(sres);
	Num_free_ica_result(ires);
	Num_free_mat(A);
}

void
test17()
{
	RecImage			*img_m, *img_p, *img, *img_slc;
	RecLoop				*slcLoop;	// cardiac phase
	RecLoop				*timeLoop;
	RecLoop				*bLoop;
	NSString			*base = @"../toshiba_images/DWI-rinkan-1024";
	int					bval = 100;		// 100, 200, 500
	NSString			*path;
	float				pscale = M_PI / 10000;
	int					i;

// mag
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-mag.img", base, bval, bval];
	img_m = [RecImage imageWithKOImage:path];
	path = [NSString stringWithFormat:@"%@/b%d-proc/b%d-phs.img", base, bval, bval];
	img_p = [RecImage imageWithKOImage:path];
	[img_p multByConst:pscale];
	[img_m makeComplexWithPhs:img_p];


//	path = [NSString stringWithFormat:@"%@/b%d-cpx.img", base, bval];
//	[img_m saveAsKOImage:path];

	bLoop  = [RecLoop loopWithDataLength:2];
	slcLoop  = [RecLoop loopWithDataLength:3];
	timeLoop = [RecLoop loopWithDataLength:128];
	img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:timeLoop, bLoop, slcLoop, [img_m yLoop], [img_m xLoop], nil];
	[img copyImageData:img_m];
	[img dumpLoops];
	path = [NSString stringWithFormat:@"%@/b%d.img", base, bval];
	[img saveAsKOImage:path];

	img = [img sliceAtIndex:1 forLoop:bLoop];
	for (i = 0; i < 3; i++) {
		img_slc = [img sliceAtIndex:i forLoop:slcLoop];

		// baseline phase removal for each slice
		img_p = [img_slc avgForLoop:timeLoop];
		[img_p toUnitImage];
		[img_p conjugate];
		[img_slc multByImage:img_p];

		// 0th phase removal for each time point
		[img_slc pcorr];
		
		// unwrap
		img_p = [img_slc unwrap2d];
		// ###

		path = [NSString stringWithFormat:@"%@/b%d-%d.img", base, bval, i];
		[img_slc saveAsKOImage:path];
		path = [NSString stringWithFormat:@"%@/b%d-%d-p.img", base, bval, i];
		[img_p saveAsKOImage:path];
	}
}

// * line gridding
// * correlation in k-space
// * start from raw only
// * shift correction
// * re-gridding
// * iterative
// * termination condition
//  return value (cumulative)

void
dump_sft(float *dx, float *dy, float *dp, int n)
{
	int		i;
	for (i = 0; i < n; i++) {
		printf("%d %f %f %f\n", i, dx[i], dy[i], dp[i]);
	}
}

void
test18()	// 2D k-traj correction
{
	RecImage	*raw, *raw_avg;
	RecLoop		*ch, *sl, *xLp, *yLp;
	RecImage	*img, *img2, *img3;
	RecImage	*tmp_img;
	RecImage	*theta;
	RecImage	*traj, *traj2, *raw2;
	RecGridder	*grid;
	RecImage	*corr;
	RecLoopControl	*lc;
	int			i, nView, j, nPt, n;
	NSPoint		pt;
	RecImage	*kxSft, *kySft, *phs;
	float		*xs, *ys, *ps;
	RecImage	*dkxSft, *dkySft, *dphs;
	float		*dx, *dy, *dp, w = 0.1;
	float		*kx, *ky;
	float		*re, *im;
	float		err, err0 = 0;
	int			iter;
	int			cropDim = 64;

// initial shift correction is necessary
//	raw = [RecImage imageFromFile:@"optical/raw1.recimg" relativePath:YES];		// before shift correction
//	raw = [RecImage imageFromFile:@"optical/raw2.recimg" relativePath:YES];		// after shift correction
	raw = [RecImage imageFromFile:@"optical/raw3.recimg" relativePath:YES];		// after 0th correction
//	raw = [RecImage imageWithKOImage:@"optical/raw_avg.img"];	// scaled down avg
	theta = [RecImage imageWithKOImage:@"optical/IMG_gasort"];

// make avg version
	ch = [RecLoop findLoop:@"Channel"];
	sl = [RecLoop findLoop:@"Slice"];
	xLp = [RecLoop loopWithDataLength:[raw xDim]];
	yLp = [RecLoop loopWithDataLength:[raw xDim]];
	raw_avg = [raw avgForLoop:ch];
	raw_avg = [raw_avg avgForLoop:sl];
	[raw_avg crop:[raw_avg xLoop] to:cropDim];

	[raw_avg saveAsKOImage:@"raw_avg.img"];

	traj  = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:[raw_avg yLoop], [raw_avg xLoop], nil];
	[traj initRadialTrajWithTab:theta startPoint:0];
	nView = [traj yDim];
	nPt = [traj xDim];

	img  = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nPt yDim:nPt];

	img2 = [RecImage imageWithImage:img];
	img3 = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[traj yLoop], [img yLoop], [img xLoop], nil];

// cumulative (return value)
	kxSft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	kySft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	phs  = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	xs = [kxSft data];
	ys = [kySft data];
	ps = [phs data];
	[kxSft clear];
	[kySft clear];
	[phs   clear];

// delta
	dkxSft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	dkySft = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	dphs  = [RecImage imageOfType:RECIMAGE_REAL withLoops:[traj yLoop], nil];
	dx = [dkxSft data];
	dy = [dkySft data];
	dp = [dphs data];

	fft_dbg = NO;
// ========================== iteration loop
	for (iter = 0; iter < 20; iter++) {
	// make ref
	//	grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:YES];
		grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:NO];
		// density correction makes iteration unstable ... apply correction in final gridding only
		[grid grid1Data:raw_avg withControl:[raw_avg control] to:img];
	//	[grid grid2d:raw to:img2];
		[img saveAsKOImage:@"IMG_ref_k"];
		tmp_img = [img copy];
		[tmp_img fft2d:REC_FORWARD];
		[tmp_img saveAsKOImage:@"IMG_ref"];

		// === line gridding
		for (i = 0; i < nView; i++) {
			// traj of one line
			traj2 = [traj sliceAtIndex:i forLoop:[traj yLoop]];
			raw2 = [raw_avg sliceAtIndex:i forLoop:[raw_avg yLoop]];
			lc = [raw2 control];
			[lc deactivateX];
			// grid one line
			grid = [RecGridder gridderWithTrajectory:traj2 andRecDim:[img xDim] nop:1 densityCorrection:NO];
			[img2 clear];
			[grid grid1Data:raw2 withControl:lc to:img2];
			// copy to multi-slice image
			[img3 copySlice:img2 atIndex:i];
		}
		[img3 saveAsKOImage:@"IMG_lines"];
		corr = [img3 xyCorrelationWith:img width:0.2 triFilt:NO];
		[corr saveAsKOImage:@"IMG_corr"];
		for (i = 0; i < nView; i++) {
			// pixels, radian
			pt = [[corr sliceAtIndex:i] findPeak2DwithPhase:dp + i];
			dx[i] = pt.x;
			dy[i] = pt.y;
		}

		[dkxSft gauss1DLP:w forLoop:[dkxSft xLoop]];
		[dkySft gauss1DLP:w forLoop:[dkySft xLoop]];
		[dphs  gauss1DLP:w forLoop:[dphs xLoop]];
		dx = [dkxSft data];
		dy = [dkySft data];
		dp = [dphs data];

		err = 0;
		for (i = 0; i < nView; i++) {
			err += dx[i]*dx[i] + dy[i]*dy[i];
			xs[i] += dx[i];
			ys[i] += dy[i];
			ps[i] += dp[i];
		}
		printf("%d %f\n", iter, err);
		if (iter == 0) {
			err0 = err + 1;
		}
		if (err >= err0) {
			break;
		}
		err0 = err;

		// === kxy, p0 correcion
		kx = [traj real];
		ky = [traj imag];
		re = [raw_avg real];
		im = [raw_avg imag];
		for (i = 0; i < nView; i++) {
			for (j = 0; j < nPt; j++) {
				kx[i * nPt + j] -= dx[i] / nPt;
				ky[i * nPt + j] -= dy[i] / nPt;
			}
			Rec_corr0th(re + i * nPt, im + i * nPt, 1, nPt, dp[i]);
		}
	//	if (iter == 0) {
	//		dump_sft(dx, dy, dp, nView);
	//	}
	}
// ========================== iteration loop

//	dump_sft(dx, dy, dp, nView);
	dump_sft(xs, ys, ps, nView);

// make corrected data (traj, raw)
	traj  = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:[raw yLoop], [raw xLoop], nil];
	[traj initRadialTrajWithTab:theta startPoint:0];
	nView = [traj yDim];
	nPt = [traj xDim];

	// === kx, ky correcion
	kx = [traj real];
	ky = [traj imag];
	for (i = 0; i < nView; i++) {
		for (j = 0; j < nPt; j++) {
			kx[i * nPt + j] -= xs[i] / nPt;
			ky[i * nPt + j] -= ys[i] / nPt;
		}
	}
	// === p0 correction #####
[raw saveAsKOImage:@"IMG_before"];
	n = [raw yDim];
	for (i = 0; i < nView; i++) {
		img2 = [raw sliceAtIndex:i forLoop:[raw yLoop]];
		re = [img2 real];
		im = [img2 imag];
		Rec_corr0th(re, im, 1, [img2 dataLength], dp[i]);
		[raw copySlice:img2 atIndex:i forLoop:[raw yLoop]];
	}
[raw saveAsKOImage:@"IMG_after"];

// final regridding (start point for motion correction)
	img  = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, yLp, xLp, nil];
	grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
	[grid grid2d:raw to:img];
	[img saveAsKOImage:@"IMG_final.img"];
	img2 = [img combineForLoop:ch];
	[img2 saveAsKOImage:@"IMG_c.img"];
}

// 117
void
test19()	// plot laplacian with shift scale
{
//	RecImage	*img = [RecImage imageWithKOImage:@"imgLapScl.img"];
	RecImage	*img = [RecImage imageWithKOImage:@"imgs_scl.sav"];	// coro
	RecImage	*img2;
	RecImage	*slc;
	RecLoop		*slcLp, *sftLp;
	int			n, nSlc, nSft = 20;

	n = [img zDim];
	nSlc = n / nSft;

//printf("nSlc = %d\n", nSlc);
//exit(0);
	slcLp = [RecLoop loopWithDataLength:nSlc];
	sftLp = [RecLoop loopWithDataLength:nSft];
	img2 = [RecImage imageOfType:RECIMAGE_REAL withLoops:slcLp, sftLp, [img yLoop], [img xLoop], nil];
	[img2 copyImageData:img];

	[img2 swapLoop:[img yLoop] withLoop:slcLp];

[img2 laplace1dForLoop:[img2 yLoop] direction:REC_FORWARD];
[img2 magnitude];
[img2 gauss2DLP:0.05];
[img2 saveAsKOImage:@"imgs_lapy.img"];
exit(0);


	[img2 swapLoop:sftLp withLoop:slcLp];
	[img2 saveAsKOImage:@"imgLapScl_reorder.img"];
	slc = [img2 sliceAtIndex:116 forLoop:slcLp];
	[slc saveAsKOImage:@"imgLapScl_slc.img"];

	img = [RecImage imageWithKOImage:@"imgs.img"];
	[img2 swapLoop:sftLp withLoop:slcLp];
	[img2 copyImageData:img];
	[img2 swapLoop:sftLp withLoop:slcLp];
	[img2 saveAsKOImage:@"imgs_reorder.img"];
}

void
test20()	// test routine for maxIndex step
{
	RecImage	*img, *img_f;
	RecImage	*lap = [RecImage imageWithKOImage:@"imgLapScl.img"];
	RecImage	*mxSft;
	RecImage	*tmp;	// tmp image without loop hierarchy -> copy image data to final destination
	RecLoop		*slcLp, *sftLp;
	int			n, nSlc, nSft = 10;

	tmp = [RecImage imageWithKOImage:@"imgs.img"];
	n = [tmp zDim];
	nSlc = n / nSft;
	slcLp = [RecLoop loopWithDataLength:nSlc];
	sftLp = [RecLoop loopWithDataLength:nSft];

	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:sftLp, slcLp, [tmp yLoop], [tmp xLoop], nil];
	[img copyImageData:tmp];

	tmp = [RecImage imageWithKOImage:@"imgLapScl.img"];
	lap = [RecImage imageWithImage:img];
	[lap copyImageData:tmp];


	mxSft = [lap peakIndexForLoop:sftLp];	// sharper with larger laplacian
	[mxSft saveAsKOImage:@"imgMxSft.img"];

// select best image
	img_f = [img selectSft:mxSft];
	[img_f saveAsKOImage:@"img_final.img"];



}

void
test21()	// coil combine with thres
{
	RecImage	*imgr = [RecImage imageWithKOImage:@"../toshiba_images/2018-02-19/15/img.img"];
	RecImage	*img;
	RecLoop		*zLp, *chLp, *xLp, *yLp;

	system("rm *.img");
	system("rm IMG_*");
	xLp = [RecLoop loopWithDataLength:[imgr xDim]];
	yLp = [RecLoop loopWithDataLength:[imgr yDim]];
	zLp = [RecLoop loopWithDataLength:72];
	chLp = [RecLoop loopWithDataLength:24];
	img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:chLp, zLp, yLp, xLp, nil];
	[img copyImageData:imgr];
	img = [img combineForLoop:chLp threshold:0.5];
	[img saveAsKOImage:@"IMG_comb.img"];
}

void
test22()
{
	RecImage		*img, *mask, *a_im, *u_im, *vt_im;
	RecImage		*prof, *comb;
	RecLoopControl	*lc;
	RecLoop			*ch, *sl;
	int				xDim, yDim, nCh;
	RecImage		*tmp_im;
	Num_svd_result	*sres;
	Num_mat			*A;
	int				nr, nc;
	float			*p1, *q1;	// img
	float			*p2, *q2;	// A
	float			*p3, *q3;	// prof
	float			*p4, *q4;	// U
	float			*p5, *q5;	// comb
	float			*p6, *q6;	// Vt
	int				i, j, k;	// position of UL corner
	int				d = 8;
	int				ii, jj;	// index into block
	NSString		*path;

// rtf3d volume
	system("rm *.img");
	system("rm IMG_*");
	img = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [img zLoop];
	img = [img sliceAtIndex:10 forLoop:sl];	// 32
//	img = [RecImage imageWithKOImage:@"/Users/oshio/Math/SMASH/spine_coil_phantom/IMG.cpx"];
//	ch = [img zLoop];

//[img fft2d:REC_INVERSE];
//[img takeEvenLines];
//[img fft2d:REC_FORWARD];
	[img saveAsKOImage:@"img.img"];

	xDim = [img xDim];
	yDim = [img yDim];
	nCh = [ch dataLength];

	nr = d * d;
	nc = nCh;

	prof = [RecImage imageOfType:RECIMAGE_COMPLEX withImage:img];
	p3 = [prof real];
	q3 = [prof imag];
	lc = [img control];
	[lc removeLoop:ch];
	comb = [RecImage imageOfType:RECIMAGE_COMPLEX withControl:lc];
	p5 = [comb real];
	q5 = [comb imag];

// SVD-based est
	p1 = [img real];
	q1 = [img imag];
	a_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nr yDim:nc];
	tmp_im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:d yDim:d zDim:nc];
	p2 = [a_im real];
	q2 = [a_im imag];
	for (i = 0; i < yDim - d; i++) {
		printf("%d\n", i);
		for (j = 0; j < xDim - d; j++) {
			// svd
			for (ii = 0; ii < d; ii++) {
				for (jj = 0; jj < d; jj++) {
					for (k = 0; k < nCh; k++) {
						p2[k * nr + ii * d + jj] = p1[k * xDim * yDim + (i + ii) * xDim + (j + jj)];
						q2[k * nr + ii * d + jj] = q1[k * xDim * yDim + (i + ii) * xDim + (j + jj)];
					}
				}
			}
			A = Num_im_to_m(a_im);
			sres = Num_svd(A);
			u_im = Num_m_to_im(sres->U);
			[u_im trans];
			p4 = [u_im real];
			q4 = [u_im imag];

			vt_im = Num_m_to_im(sres->Vt);
			p6 = [vt_im real];
			q6 = [vt_im imag];
	//	[vt_im saveAsKOImage:@"IMG_vt"];
	//		path = [NSString stringWithFormat:@"U-%d-%d.img", i, j];
	//		[u_im saveAsKOImage:path];
	//		[tmp_im copyImageData:vt_im];
	//		path = [NSString stringWithFormat:@"E-%d-%d.img", i, j];
	//		[tmp_im saveAsKOImage:path];

			// copy to profile image
			int s = 0;
			for (k = 0; k < nCh; k++) {
				p3[k * xDim * yDim + (i + ii - d/2) * xDim + (j + jj - d/2)] = p4[s * nc + k];
				q3[k * xDim * yDim + (i + ii - d/2) * xDim + (j + jj - d/2)] = q4[s * nc + k];
			}
			p5[i * xDim + j] = p6[d/2 * d + d/2];
			q5[i * xDim + j] = q6[d/2 * d + d/2];
		}
		Num_free_mat(A);
		Num_free_svd_result(sres);
	}
	[prof saveAsKOImage:@"IMG_prof"];
	[comb saveAsKOImage:@"IMG_comb"];
//	[prof gauss2DLP:0.1];
//	[prof saveAsKOImage:@"IMG_prof_f"];

//	comb = [prof combineForLoop:ch];
//	[comb saveAsKOImage:@"IMG_coil_comb"];

	[prof conjugate];
	[img multByImage:prof];
	[img saveAsKOImage:@"IMG_pcorr"];
	img = [img sumForLoop:ch];
	[img saveAsKOImage:@"IMG_ccomb"];
	

}

// im1 -> (im1 . im2) / |im1||im2| (vector in channel direction)
// vDSP version slower -> NSOperation ?
void
err_proj(RecImage *proj, RecImage *im1, RecImage *im2, RecLoop *lp)
{
	RecLoopControl	*lc;
	int				i, j, n, len;
	int				ix, skip;
	float			*p0, *q0;
	float			*p1, *p2;
	float			*q1, *q2;
	float			mg1, mg2, pr, pi;
//	DSPSplitComplex	A, B, C;
//	float			ssq;

	lc = [im1 control];
	[lc deactivateLoop:lp];
	n = [lc loopLength];
	len = [lp dataLength];
	skip = [im1 skipSizeForLoop:lp];
	p0 = [proj real];
	q0 = [proj imag];
	for (i = 0; i < n; i++) {
		p1 = [im1 currentDataWithControl:lc];
		q1 = p1 + [im1 dataLength];
		p2 = [im2 currentDataWithControl:lc];
		q2 = p2 + [im2 dataLength];

//		A.realp = p1;	A.imagp = q1;
//		B.realp = p2;	B.imagp = q2;
//		C.realp = &pr;	C.imagp = &pi;
//		vDSP_zdotpr(&A, skip, &B, skip, &C, len);
//		mg1 = 0;
//		vDSP_svesq(p1, skip, &ssq, len); mg1 += ssq;
//		vDSP_svesq(q1, skip, &ssq, len); mg1 += ssq;
//		mg2 = 0;
//		vDSP_svesq(p2, skip, &ssq, len); mg2 += ssq;
//		vDSP_svesq(q2, skip, &ssq, len); mg2 += ssq;

		mg1 = mg2 = 0;
		pr = pi = 0;
		for (j = ix = 0; j < len; j++, ix += skip) {
			mg1 += p1[ix]*p1[ix] + q1[ix]*q1[ix];
			mg2 += p2[ix]*p2[ix] + q2[ix]*q2[ix];
			pr  += p1[ix]*p2[ix] + q1[ix]*q2[ix];
			pi  += p2[ix]*q1[ix] - p1[ix]*q2[ix];
		}
		mg1 = sqrt(mg1);
		mg2 = sqrt(mg2);
	//	p0[i] = pr/mg1;
	//	q0[i] = pi/mg1;
		p0[i] = pr/mg1/mg2;
		q0[i] = pi/mg1/mg2;
		[lc increment];
	}
}

RecImage *
mag_mask(RecImage *m, float thres)
{
	int		i, n;
	float	*p, *q;
	float	mx;
	RecImage	*msk = [m copy];

	mx = [m maxVal];
	thres *= mx;
	n = [msk dataLength];
	p = [m data];
	q = [msk data];

	for (i = 0; i < n; i++) {
		if (p[i] > thres) {
			q[i] = 1.0;
		} else {
			q[i] = p[i] / thres;
		}
	}
	return msk;
}

// updated (temporary) version
// make this consitent with doc
// 11-18-2019
void
test23()
{
	RecImage		*y, *yy, *s, *ss, *m, *mm, *err;
	RecImage		*mask;
	RecLoopControl	*lc;
	RecLoop			*ch, *sl;
	int				nCh, n;
	int				iter;
	NSString		*path;
	float			w = 8.0 / 256; // 0.03
//	float			gain = 0.5;
	float			nz = 0.04;		// frac of max val
	float			en, ep, ee;
	float			tol = 1e-10;
	BOOL			threeD = NO;

// rtf3d volume
	system("rm *.img");
	system("rm IMG_*");
//	y = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	y = [RecImage imageFromFile:@"../toshiba_images/IMGsav2.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];
	sl = [y zLoop];
	if (!threeD) {
		y = [y sliceAtIndex:42 forLoop:sl];	// 32
	}

	[y saveAsKOImage:@"img.img"];

	nCh = [ch dataLength];
	n = [y dataLength];

	lc = [y control];
	[lc removeLoop:ch];

	// == (0) initial est (S-o-S)
	m = [y combineForLoop:ch];
	[m saveAsKOImage:@"IMG_m_0"];
	mask = mag_mask(m, nz);
//	mask = [m copy];
//	[mask thresAt:nz];
	[mask saveAsKOImage:@"IMG_msk"];
	[mask gauss2DLP:w];
	[mask saveAsKOImage:@"IMG_msk_f"];
//exit(0);
	ss = [y copy];
	[ss divImage:m];	// raw profile
//	[ss multByImage:mask];
	[ss saveAsKOImage:@"IMG_prof_hi"];

	ep = [y rmsVal];
	for (iter = 1; iter < 20; iter++) {
		s = [ss copy];
		// == (1) LPF
		if (threeD) {
			[s gauss3DLP:w];
		//	[s rect3DLP:w];
		} else {
			[s gauss2DLP:w];
		//	[s rect2DLP:w];
		}
		// == (2) normalize
		[s normalize2ForLoop:ch withMask:mask];
//		[s normalize2ForLoop:ch];
		path = [NSString stringWithFormat:@"IMG_prof%d", iter];
		[s saveAsKOImage:path];

		// == (3) update m
		ss = [s copy];
		[ss conjugate];
		yy = [y copy];
		[yy multByImage:ss];
		m = [yy sumForLoop:ch];
		path = [NSString stringWithFormat:@"IMG_m%d", iter];
		[m saveAsKOImage:path];

		// === save coil est
		yy = [s copy];
		[yy multByImage:m];
		path = [NSString stringWithFormat:@"IMG_coil%d", iter];
		[yy saveAsKOImage:path];

		// === error calc
		err = [yy copy];
		[err subImage:y];
		path = [NSString stringWithFormat:@"IMG_diff%d", iter];
		[err saveAsKOImage:path];
		en = [err rmsVal];
		ee = (ep - en) / ep;
	//	printf("%d %f %f\n", iter, en, ee);
		printf("%d %f\n", iter, en);
		if (ee < tol) break;
		ep = en;

		// === (4) update si
		ss = [y copy];
		mm = [m copy];
		[ss multByImage:mask];
		[ss cpxDivImage:mm];
		path = [NSString stringWithFormat:@"IMG_prof_hi%d", iter];
		[ss saveAsKOImage:path];
	}

	[m saveAsKOImage:@"IMG_final"];
	[s saveAsKOImage:@"IMG_prof_final"];
}

// #### make lib method when done
// 
RecImage *
fit_with_mask(RecImage *img, RecImage *msk, int w)
{
	RecImage	*x0, *y, *dy;
	int			i;
	float		gain = 2.0;
	BOOL		threeD = ([msk dim] > 2);

	x0 = [img copy];	// initial
	[x0 multByImage:msk];
	y = [x0 copy];
	if (threeD) {
		[y cos3DLPc:w];
	} else {
		[y cos2DLPc:w];
	}
	for (i = 0; i < 2; i++) {
		dy = [y copy];
		[dy subImage:x0];
		[dy multByImage:msk];
		if (threeD) {
			[dy cos3DLPc:w];
		} else {
			[dy cos2DLPc:w];
		}
		[dy multByConst:gain];
		[y subImage:dy];
	}

	return y;
}

// new version (flexible data size, either full or part of k-space)
// ## define API and move to lib
//
//	
//
void
test24()
{
	RecImage		*x, *y, *yy, *s, *ss, *m, *err;
	RecImage		*mask;
	RecLoopControl	*lc;
	RecLoop			*ch, *sl;
	int				nCh, n;
	int				iter;
	NSString		*path;
	int				xDim0, zDim0;
	int				dx = 16;
	int				w = 4;
	float			nz = 0.05;		// 0.04 of max val
	float			en, ep, ee;
	float			tol = 1e-3;
	BOOL			threeD = NO;

// rtf3d volume
	system("rm *.img");
	system("rm IMG_*");
//	x = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
	x = [RecImage imageFromFile:@"../toshiba_images/IMGsav2.recimg" relativePath:YES];
	ch = [RecLoop findLoop:@"Channel"];

	sl = [x zLoop];
	if (threeD) {
		xDim0 = [x xDim];
		zDim0 = [x zDim];
	} else {
		x = [x sliceAtIndex:42 forLoop:sl];	// 32
		xDim0 = [x xDim];
		zDim0 = 1;
	}
	[x saveAsKOImage:@"img.img"];

// crop x -> y (simulate central  k-space acq)
	if (threeD) {
		y = [x copy];
		[y fft3d:REC_INVERSE];
		[y crop:[y xLoop] to:dx];
		[y crop:[y yLoop] to:dx];
		[y crop:[y zLoop] to:dx];
		[y fft3d:REC_FORWARD];
	} else {
		y = [x copy];
		[y fft2d:REC_INVERSE];
		[y crop:[y xLoop] to:dx];
		[y crop:[y yLoop] to:dx];
		[y fft2d:REC_FORWARD];
	}
	[y saveAsKOImage:@"img_center.img"];

	nCh = [ch dataLength];
	n = [y dataLength];

	lc = [y control];
	[lc removeLoop:ch];

// == (0) initial est (S-o-S)
	m = [y combineForLoop:ch];
	[m saveAsKOImage:@"IMG_m_0"];

// == (0) make mask
	mask = [m copy];
	[mask thresAt:nz];
	[mask saveAsKOImage:@"IMG_msk"];
	[mask gauss2DLP:(float)w/dx];
	[mask limitToVal:0.5];
	[mask saveAsKOImage:@"IMG_MSK_f"];

	// initial error
	ep = [y rmsVal];

	[m makeComplex];

	for (iter = 1; iter < 50; iter++) {
		// === (1) update si
		s = [y copy];
		[m multByImage:mask];
		[s cpxDivImage:m];
//		[s limitToVal:1.0];

		path = [NSString stringWithFormat:@"IMG_prof_hi%d", iter];
		[s saveAsKOImage:path];

		// == (2) LPF
		ss = fit_with_mask(s, mask, w);

		// == (3) normalize
		[ss normalize2ForLoop:ch];
//		[ss normalize2ForLoop:ch withMask:mask];
		path = [NSString stringWithFormat:@"IMG_prof%d", iter];
		[ss saveAsKOImage:path];

		// == (4) update m
		s = [ss copy];
		[s conjugate];
		yy = [y copy];
		[yy multByImage:s];
		m = [yy sumForLoop:ch];

		// == (5) update m
		if (threeD) {
			[m cos3DLPp:w];
		} else {
			[m cos2DLPp:w];
		}
		path = [NSString stringWithFormat:@"IMG_m%d", iter];
		[m saveAsKOImage:path];

		// === save coil est
		yy = [ss copy];
		[yy multByImage:m];
		path = [NSString stringWithFormat:@"IMG_coil%d", iter];
		[yy saveAsKOImage:path];

		// === (6) error calc
		err = [yy copy];
		[err subImage:y];
		path = [NSString stringWithFormat:@"IMG_diff%d", iter];
		[err saveAsKOImage:path];
		en = [err rmsVal];
		ee = (ep - en) / ep;
	//	printf("%d %f %f\n", iter, en, ee);
		printf("%d %f\n", iter, en);
		if (ee < tol) break;
		ep = en;

	}

// (7) zerofil to original dim
	if (threeD) {
		[ss fft3d:REC_INVERSE];
		[ss zeroFill:[ss xLoop] to:xDim0];
		[ss zeroFill:[ss yLoop] to:xDim0];
		[ss zeroFill:[ss zLoop] to:zDim0];
		[ss fft3d:REC_FORWARD];

		[mask fft3d:REC_INVERSE];
		[mask zeroFill:[mask xLoop] to:xDim0];
		[mask zeroFill:[mask yLoop] to:xDim0];
		[mask zeroFill:[mask zLoop] to:zDim0];
		[mask fft3d:REC_FORWARD];
	} else {
		[ss fft2d:REC_INVERSE];
		[ss zeroFill:[ss xLoop] to:xDim0];
		[ss zeroFill:[ss yLoop] to:xDim0];
		[ss fft2d:REC_FORWARD];

		[mask fft2d:REC_INVERSE];
		[mask zeroFill:[mask xLoop] to:xDim0];
		[mask zeroFill:[mask yLoop] to:xDim0];
		[mask fft2d:REC_FORWARD];
	}

// (8) normalize
//	[ss normalize2ForLoop:ch];
	[ss normalize2ForLoop:ch withMask:mask];
	[ss saveAsKOImage:@"IMG_prof_final"];

// (9) update m (final)
	[ss conjugate];
	yy = [x copy];
	[yy copyLoopsOf:ss];
	[yy multByImage:ss];
	m = [yy sumForLoop:ch];
	[m saveAsKOImage:@"IMG_final"];
}

// Q&D check of radial "parallel imaging"
// - make profile estimation method
// - make numerical phantom
void
test25()
{
	RecImage		*x, *y, *yy, *s, *ss, *m, *err;
	RecImage		*ref, *mg1, *mg2, *tmp;
	RecImage		*mask;
	RecLoopControl	*lc;
	RecLoop			*ch, *sl, *xLp, *yLp;
	int				nCh, n, nSlc = 72;
	int				iter;
	NSString		*path;
	int				xDim0, zDim0;
	int				dx = 64;
	int				w = 4;
	float			nz = 0.06;		// 0.04 of max val
	float			en, ep, ee;
	float			tol = 1e-3;
	BOOL			threeD = NO;

// rtf3d volume
	system("rm *.img");
	system("rm IMG_*");
//	x = [RecImage imageFromFile:@"../toshiba_images/IMGsav.recimg" relativePath:YES];
//	x = [RecImage imageFromFile:@"../toshiba_images/IMGsav2.recimg" relativePath:YES];
	ref = [RecImage imageWithKOImage:@"single_sav/img_488.img"];
//	x = [RecImage imageWithKOImage:@"single_sav/img_244.img"];
	x = [RecImage imageWithKOImage:@"single_sav/img_162.img"];

// ### make loops
	nCh = [x zDim];
	nCh /= nSlc;
	xLp = [x xLoop];
	yLp = [x yLoop];
	ch = [RecLoop loopWithDataLength:nCh];
	sl = [RecLoop loopWithDataLength:nSlc];
	tmp = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, yLp, xLp, nil];
	[tmp copyImageData:x];
	x = [tmp copy];
	[tmp copyImageData:ref];
	ref = [tmp copy];
	
//[x subImage:ref];
//[x saveAsKOImage:@"IMG_sub"];
//exit(0);

	[x copyLoopsOf:ref];
	ch = [RecLoop findLoop:@"Channel"];
	ch = [ref topLoop];

	sl = [x zLoop];
	if (threeD) {
		xDim0 = [x xDim];
		zDim0 = [x zDim];
	} else {
		x = [x sliceAtIndex:42 forLoop:sl];	// 32
		ref = [ref sliceAtIndex:42 forLoop:sl];
		xDim0 = [x xDim];
		zDim0 = 1;
	}
	[x saveAsKOImage:@"img.img"];
	tmp = [x combineForLoop:ch];
	[tmp saveAsKOImage:@"IMG_sos"];

// crop x -> y (simulate central  k-space acq)
	if (threeD) {
		y = [x copy];
		[y fft3d:REC_INVERSE];
		[y crop:[y xLoop] to:dx];
		[y crop:[y yLoop] to:dx];
		[y crop:[y zLoop] to:dx];
		[y fft3d:REC_FORWARD];
	} else {
		y = [x copy];
		[y fft2d:REC_INVERSE];
		[y crop:[y xLoop] to:dx];
		[y crop:[y yLoop] to:dx];
		[y fft2d:REC_FORWARD];
	}
	[y saveAsKOImage:@"img_center.img"];

	nCh = [ch dataLength];
	n = [y dataLength];

	lc = [y control];
	[lc removeLoop:ch];

// == (0) initial est (S-o-S)
	m = [y combineForLoop:ch];
	[m saveAsKOImage:@"IMG_m_0"];
	mask = [m copy];
	[mask thresAt:nz];
	[mask saveAsKOImage:@"IMG_msk"];
	[m makeComplex];

	// initial error
	ep = [y rmsVal];

	for (iter = 1; iter < 50; iter++) {
		// === (1) update si
		s = [y copy];
		[m multByImage:mask];
		[s cpxDivImage:m];
//		[s limitToVal:1.0];

		path = [NSString stringWithFormat:@"IMG_prof_hi%d", iter];
		[s saveAsKOImage:path];

		// == (2) LPF
		ss = fit_with_mask(s, mask, w);

		// == (3) normalize
		[ss normalize2ForLoop:ch];
		path = [NSString stringWithFormat:@"IMG_prof%d", iter];
		[ss saveAsKOImage:path];

		// == (4) update m
		s = [ss copy];
		[s conjugate];
		yy = [y copy];
		[yy multByImage:s];
		m = [yy sumForLoop:ch];

		// == (5) update m
		if (threeD) {
			[m cos3DLPp:w];
		} else {
			[m cos2DLPp:w];
		}
		path = [NSString stringWithFormat:@"IMG_m%d", iter];
		[m saveAsKOImage:path];

		// === save coil est
		yy = [ss copy];
		[yy multByImage:m];
		path = [NSString stringWithFormat:@"IMG_coil%d", iter];
		[yy saveAsKOImage:path];

		// === (6) error calc
		err = [yy copy];
		[err subImage:y];
		path = [NSString stringWithFormat:@"IMG_diff%d", iter];
		[err saveAsKOImage:path];
		en = [err rmsVal];
		ee = (ep - en) / ep;
	//	printf("%d %f %f\n", iter, en, ee);
		printf("%d %f\n", iter, en);
		if (ee < tol) break;
		ep = en;

	}

// (7) zerofil to original dim
	if (threeD) {
		[ss fft3d:REC_INVERSE];
		[ss zeroFill:[ss xLoop] to:xDim0];
		[ss zeroFill:[ss yLoop] to:xDim0];
		[ss zeroFill:[ss zLoop] to:zDim0];
		[ss fft3d:REC_FORWARD];
	} else {
		[ss fft2d:REC_INVERSE];
		[ss zeroFill:[ss xLoop] to:xDim0];
		[ss zeroFill:[ss yLoop] to:xDim0];
		[ss fft2d:REC_FORWARD];
	}

// (8) normalize
	[ss normalize2ForLoop:ch];
	[ss saveAsKOImage:@"IMG_prof_final"];

// (9) update m (final)
	[ss conjugate];
	yy = [x copy];
	[yy copyLoopsOf:ss];
	[yy multByImage:ss];
	m = [yy sumForLoop:ch];
	[m saveAsKOImage:@"IMG_final"];


// dif
	tmp = [ref copy];
	[tmp combineForLoop:ch];
	[m combineForLoop:ch];
	[tmp copyLoopsOf:m];
	[m subImage:tmp];
	[m saveAsKOImage:@"IMG_dif"];
	printf("RMS = %f\n", [m rmsVal]);
}

void
test_sct(RecImage *tgt, RecImage *ref)
{
	float	*p, *q;
	int		i, n = [tgt xDim];
	p = [tgt data];
	q = [ref data];
	for (i = 0; i < n; i++) {
		printf("%f %f\n", q[i], p[i]);
	}
}

// shift est using central view
// phantom OK
// lib method ok - (RecImage *)shiftFromK0;	// shift estimation from center-of-kspace view (POCS)

void
test26()
{
	RecImage	*prj = [RecImage imageWithKOImage:@"central_nav_t.img"];	// time order (now multi-channel)
	RecImage	*nav = [RecImage imageWithKOImage:@"nav_sft_t.img"];
//	RecImage	*prj = [RecImage imageWithKOImage:@"central_nav.img"];		// angle order (now multi-channel)
//	RecImage	*nav = [RecImage imageWithKOImage:@"nav_sft.img"];
	RecImage	*mean;
	RecImage	*st, *mv, *ms, *sft, *nsft;		// stationary[xDim], moving[xDim], shift[yDim]
//	RecImage	*win;
//	int			wid = 4;
	float		nz = 0.005;	// 0.001
	RecImage	*est, *dif;
	RecLoop		*tLp, *xLp, *chLp;
	int			i, j, k, iter, nIter;
	int			xDim, yDim, nCh;
	float		*m, *p, *s, *sf;
	float		x, val, err, mx;
	NSString	*path;

	printf("test 26\n");
//	system("rm IMG_*.img");

	[nav negate];
	[nav saveAsKOImage:@"IMG_ref.img"];


// method test
if (1) {
	prj = [prj combineForLoop:[prj zLoop]];
	sft = [prj shiftFromK0];
	test_sct(sft, nav);
	exit(0);
}

	[prj removePointLoops];
	[prj takeRealPart];

	if (0) {	// make test data
		int tDim = 256;
		int	xDim = 64;
		prj = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim yDim:tDim];
		mv = [prj sliceAtIndex:0 forLoop:[prj yLoop]];
		ms = [mv copy];
		st = [mv copy];
		sft = [RecImage imageOfType:RECIMAGE_REAL xDim:tDim];

		s = [st data];
		m = [mv data];
		p = [prj data];
		sf = [sft data];

		p = [prj data];
		for (i = 0; i < xDim; i++) {
			x = ((float)i - xDim/2) / xDim;
			// st
			val = 0.4 - x * x;
			if (val > 0) {
				val = sqrt(val);
			} else {
				val = 0;
			}
			s[i] = val;
			// mv
			val = 0.04 - x * x;
			if (val > 0) {
				val = sqrt(val);
			} else {
				val = 0;
			}
			m[i] = val;
		}
//		[st saveAsKOImage:@"IMG_st"];
//		[mv saveAsKOImage:@"IMG_mv"];
		
		// sft
		for (i = 0; i < tDim; i++) {
			x = (float)i * M_PI * 10 / tDim;
			sf[i] = sin(x);
		}
//		[sft saveAsKOImage:@"IMG_sft"];

		for (i = 0; i < tDim; i++) {
			ms = [mv shift1dLoop:[ms xLoop] by:sf[i] * 5.0];	// +- 5 pixel shift
			m = [ms data];
			for (j = 0; j < xDim; j++) {
				val = 0.2 * m[j];
				p[i * xDim + j] = val;
			}
		}
		[prj saveAsKOImage:@"IMG_ms"];

		[prj addImage:st];
		[prj addGWN:nz relative:YES]; // 0.001, 0.01
//	[prj zeroFill:[prj xLoop] to:128];

		[prj saveAsKOImage:@"IMG_phantom"];
	}

// ============ main proc start ======
	[prj dumpLoops];
	xLp = [prj xLoop];
	tLp = [prj yLoop];
	chLp = [prj zLoop];
	xDim = [prj xDim];
	yDim = [prj yDim];
	nCh = [prj zDim];

// normalize (mag)
	mean = [prj avgForLoop:[prj xLoop]];
	m = [mean data];
	p = [prj data];
	for (k = 0; k < nCh; k++) {
		for (i = 0; i < yDim; i++) {
		//	printf("%d %f\n", i, m[i]);
			for (j = 0; j < xDim; j++) {
				if (m[i] != 0) {
					p[((k * yDim) + i) * xDim + j] /= m[i];
				}
			}
		}
	}
	[prj saveAsKOImage:@"IMG_in.img"];

	// gradient along x
	if (1) {
		int ix;
		dif = [RecImage imageWithImage:prj];
		p = [prj data];
		m = [dif data];
		for (k = 0; k < nCh; k++) {
			for (i = 0; i < yDim; i++) {
				for (j = 1; j < xDim; j++) {
					ix = ((k * yDim) + i) * xDim + j;
					m[ix] = p[ix] - p[ix - 1];
				}
			}
		}
		[dif saveAsKOImage:@"IMG_dif.img"];
		prj = [dif copy];
	}
	[prj saveAsKOImage:@"IMG_in-st.img"];
	mv = [prj sdForLoop:[prj yLoop]];
	[mv saveAsKOImage:@"IMG_sd.img"];

// tri win
	[prj fTriWin1DforLoop:[prj xLoop]];

	[prj saveAsKOImage:@"IMG_win.img"];

	mv = [prj sliceAtIndex:0 forLoop:[prj yLoop]];
	st = [RecImage imageWithImage:mv];

	[st saveAsKOImage:@"IMG_st0.img"];
	[mv saveAsKOImage:@"IMG_mv0.img"];
//exit(0);

// === iteration ===
	nIter = 20;	// 20
	for (iter = 1; iter <= nIter; iter++) {
		est = [prj copy];
		[est subImage:st];
		sft = [est corr1dWithRef:mv];	// unit:pixels
	//	sft = [[est sliceAtIndex:18 forLoop:[est zLoop]] corr1dWithRef:[mv sliceAtIndex:18 forLoop:[mv yLoop]]];

		mx = [sft meanVal];
		[sft addConst:-mx];

		nsft = [sft copy];
		[nsft negate];
		
		[est clear];
		[est copyImage:mv];
		est = [est correctShift1d:nsft forLoop:[mv xLoop]]; // pixels

		[est addImage:st];
		dif = [prj copy];
		[dif subImage:est];
		[dif fTriWin1DforLoop:[dif xLoop]];


		err = [dif rmsVal];
		printf("%d %e\n", iter, err);
		mean = [dif avgForLoop:[dif yLoop]];

	//	[mean multByConst:1.0];
		[st addImage:mean];

		mv = [prj copy];
	 	[mv subImage:st];
		mv = [mv correctShift1d:sft forLoop:[mv xLoop]];
		mv = [mv avgForLoop:[mv yLoop]];

		if (iter == nIter) {
			int ix = 1;
			test_sct(sft, nav);
			path = [NSString stringWithFormat:@"IMG_sft%d.img", ix];
			[sft saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_mv%d.img", ix];
			[mv saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_st%d.img", ix];
			[st saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_dif%d.img", ix];
			[dif saveAsKOImage:path];
			path = [NSString stringWithFormat:@"IMG_ms%d.img", ix];
			[est saveAsKOImage:path];
		}
	}
}

void
test27()	// selecting best focused shift
{
	RecImage	*img;
	RecImage	*imgs, *tmp_img, *mxSft;
	RecLoop		*xLp, *yLp, *zLp, *scLp;
//	int			i, nImg;
	int			nScl = 10;
	int			bestScl;
	float		lap_w = 0.2;

	fft_dbg = NO;

//	 [tmp_img saveAsKOImage:@"imgs_scl.img"];
	img = [RecImage imageWithKOImage:@"imgs_scl.sav"];	// y, scl, z, x
	bestScl = 4;	// 0.6

	[img dumpLoops];
	xLp   = [RecLoop loopWithDataLength:256];
	yLp   = [RecLoop loopWithDataLength:256];
	zLp   = [RecLoop loopWithDataLength:128];
	scLp = [RecLoop loopWithDataLength:nScl];
	imgs =  [RecImage imageOfType:RECIMAGE_REAL withLoops:yLp, scLp, zLp, xLp, nil];
	[imgs copyImageData:img];	// preserve
	
// (5) === find best focus scale ====
	printf("5: select best shift\n");
	tmp_img = [imgs copy];
//	[tmp_img laplace2d:REC_FORWARD];
	[tmp_img grad1dForLoop:[tmp_img yLoop]];
	[tmp_img magnitude];
	[tmp_img square];
	[tmp_img gauss3DLP:lap_w];
	[tmp_img swapLoop:scLp withLoop:[img yLoop]];	// ## chk
	[tmp_img saveAsKOImage:@"imgs_scl_lap.img"];

	mxSft = [tmp_img peakIndexForLoop:scLp];	// sharper with larger laplacian
	[mxSft saveAsKOImage:@"imgMxSft.img"];		// ok
	[imgs swapLoop:scLp withLoop:[img yLoop]];	// ## chk
	tmp_img = [imgs selectSft:mxSft];			// ### loop ? [scl z x]
	// ####
	[tmp_img saveAsKOImage:@"imgs_final.img"];

	[mxSft gauss2DLP:lap_w];
	[mxSft saveAsKOImage:@"imgMxSft_f.img"];
	tmp_img = [imgs selectSftF:mxSft];
	[tmp_img saveAsKOImage:@"imgs_final_F.img"];


	img = [tmp_img copy];
	img = [img scale1dLoop:zLp to:122];
	[img saveAsKOImage:@"imgs_scl_cor.img"];
	tmp_img = [img copy];
	[tmp_img swapLoop:xLp withLoop:yLp];
	[tmp_img saveAsKOImage:@"imgs_scl_sag.img"];
	tmp_img = [img copy];
	[tmp_img swapLoop:[tmp_img yLoop] withLoop:yLp];
	[tmp_img saveAsKOImage:@"imgs_scl_ax.img"];

	img = [imgs sliceAtIndex:0 forLoop:scLp];
	img = [img scale1dLoop:zLp to:122];
	[img saveAsKOImage:@"imgs_0_cor.img"];
	tmp_img = [img copy];
	[tmp_img swapLoop:xLp withLoop:yLp];
	[tmp_img saveAsKOImage:@"imgs_0_sag.img"];
	tmp_img = [img copy];
	[tmp_img swapLoop:[tmp_img yLoop] withLoop:yLp];
	[tmp_img saveAsKOImage:@"imgs_0_ax.img"];
}

