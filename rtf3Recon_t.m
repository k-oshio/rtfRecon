//
//	rtf3 recon (Toshiba version)
//
//	K. Oshio	1-27-2015
//
//  merge with rtf3Recon when finished
//
// === plans ===
//	* golden angle radial, M_PI * (3 - sqrt(5))
//	* slice half fourier
//	* compare with zero-fill
//	* golden angle + ref (golden angle corrected to 1 + sqrt(5))
//	* rad pcorr ok (ref based)
//	* chk echo
//	* peak detection with stationary part
//	* sort view with proj angle
//	* self rad phase corr
//	* mapless coil intensity/phase correction
//	X arm removal -> just skip lateral channel
//	X rigid correction
//	- non-rigid correction
//		*3D def est using binned reconstruction
//			rigidCorr() version has better SNR -> chk and copy
//		multi-channel, multi-position sft
//		

#import <RecKit/RecKit.h>
#import <sys/time.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

RecImage *	rigidCorr(RecImage *pwi, RecImage *img, int maxIter);
//RecImage *	nonRigidCorr(RecImage *pwi, RecImage *img, int maxIter);
void		test_diff(RecImage *pw, RecImage *img);

// global
RecLoop		*ch, *sl, *pe, *rd;
RecLoop		*kx, *ky;
RecImage	*radMap, *traj;
int			coil = Coil_None;	// TOSHIBA_15;
RecGridder	*grid;
float		thr2RR = -1;	// 0.25;
float		thrDB  = 0;		//0.5; // 0.5;	// pixels;

int
main(int ac, char *av[])
{
	NSString		*path;
	RecImage		*raw, *img, *pw, *pw_c, *pws, *pwr, *img_c;
	RecImage		*pwo, *pwi, *imgo, *imgi;
	RecImage		*gaTab = nil;
	RecImage		*gaSort;
	RecImage		*corr, *sft;
	RecImage		*tmp_img;

    RecImage        *img_coro;
	RecImage		*img_lap;	// laplacian

	int				recDim = 512;
	// recon params
	BOOL			halfFT		= NO;
	BOOL			goldenAngle	= YES;
	BOOL			sort		= YES;
	int				nRef		= 8;
	int				echoStart	= 0; //112;
	int				nop			= 1; //2;

	BOOL			rigid		= NO;
	BOOL			nonRigid	= NO;
	BOOL			hires		= NO;

TIMER_ST
    @autoreleasepool {
        system("rm *.img");
        system("rm IMG_*");
        system("rm sft*.txt");

	//	[RecImage setFFTdbg:NO];	// turn off FFT unit check

// (1) === load raw data
		if (ac < 2) {
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2015-10-22/Run272.7683.07"];	// volunteer1
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2015-10-23/Run275.7683.06"];	// volunteer2
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2015-10-23/Run275.7683.11"];	// volunteer3
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2015-11-12/Run331.7683.08"];	// volunteer
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2015-11-12/Run331.7683.09"];	// volunteer + SAT
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-02-04/Run492.7683.12"];	// half fourier (2mm)
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-02-04/Run492.7683.11"];	// short TE (3mm)
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-02-04/Run493.7683.14"];	// half echo
		// ========= new contract ======
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-12-10/Run54033.5212.11-Radial180"];	// 180
		//	goldenAngle = NO;
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2016-12-10/Run54033.5212.10-GoldenAngle"];	// g.a.
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-04-07/Run2109.6953.06-PW_RMC_OFF"];	// g.a.
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-04-07/Run2109.6953.08-PW_RMC_ON"];	// g.a.
		//	echoStart = 222; recDim = 512; nop = 1;
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-09-28/Run23531.-5213.03" traj:&gaTab];	// FFE,  1, 729, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-09-28/Run23531.-5213.04" traj:&gaTab];	// SSFP, 1, 729, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-09-28/Run23531.-5213.05" traj:&gaTab];	// FFE,  1, 600, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-09-28/Run23531.-5213.06" traj:&gaTab];	// SSFP, 1, 600, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-09-28/Run23531.-5213.07" traj:&gaTab];	// SSFP, 0.5, 600, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-10-11/Run23576.-5213.07" traj:&gaTab];	// SSFP, 0.5, 600, 1, 8
		//	raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-10-24/Run23660.-5213.08" traj:&gaTab];	// SSFP, 0.5, 600, 1, 8
			raw = [RecImage imageWithToshibaFile:@"../toshiba_images/2017-10-24/Run23660.-5213.10" traj:&gaTab];	// SSFP, 0.5, 600, 1, 8
			echoStart = 0;
		} else {
			path = [NSString stringWithFormat:@"%s", av[1]];
			raw = [RecImage imageWithToshibaFile:path];	// [ch, *sl, pe, *rd]
		}
		if (raw == nil) {
			printf("file not found [%s]\n", [path UTF8String]);
			exit(0);
		}
		[raw saveAsKOImage:@"raw0.img"];
	// raw data loops
		rd = [RecLoop findLoop:@"Read"];
		pe = [RecLoop findLoop:@"Phase"];
		ch = [RecLoop findLoop:@"Channel"];
		sl = [RecLoop findLoop:@"Slice"];
		kx = [RecLoop loopWithName:@"kx" dataLength:recDim];
		ky = [RecLoop loopWithName:@"ky" dataLength:recDim];
TIMER_END("load raw");

		sl = [raw cropToPo2:sl];		// implement dft (for small len) later ###
		if (halfFT) {
			sl = [raw halfFT2:sl];				// [ch, pe,  sl,  rd] (pw)
		} else {
			sl = [raw sftCenter:sl];
		}
//		if ([raw xDim] > 256) {
//		    [raw freqCrop];
//			rd = [raw xLoop];
//		}
	//	[raw saveAsKOImage:@"raw.img"];
	//	[raw dumpLoops];

	// remove lateral channel
	//	ch = [raw removeSliceAtIndex:14 forLoop:ch];
	//	ch = [raw removeSliceAtIndex:13 forLoop:ch];
	//	ch = [raw removeSliceAtIndex:10 forLoop:ch];
	//	ch = [raw removeSliceAtIndex: 5 forLoop:ch];
	
		// remove empty channel
		int i;
		for (i = 14; i >=4; i--) {
			ch = [raw removeSliceAtIndex:i forLoop:ch];
		}
		[raw dumpLoops];
		[raw saveAsKOImage:@"raw1.img"];
	
// (2) === pre-processing
		// create radial trajectory
		traj  = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rd, nil];
		radMap = [RecImage imageOfType:RECIMAGE_MAP  withLoops:pe, rd, nil];

		if (goldenAngle) {
			if (!gaTab) {
				gaTab = Rec_goldenAngleTab(pe, nRef);
			}
			[gaTab saveAsKOImage:@"IMG_gatab"];

			if (sort) {
				gaSort = Rec_goldenAngleSort(gaTab);
			}
			
			// reorder both raw & traj using gaSort tab
			if (sort) {	// chk
				[raw reorder:pe withTab:gaSort];
				tmp_img = [raw avgForLoop:ch];	// ### takes too long...
			//	tmp_img = [raw copy];
				[tmp_img fft1d:sl direction:REC_INVERSE];
				[tmp_img saveAsKOImage:@"raw_reorder.img"];
				[tmp_img fft1d:rd direction:REC_FORWARD];
				[tmp_img saveAsKOImage:@"sin_reorder.img"];
				[traj initRadialTrajWithSortTab:gaSort startPoint:echoStart];	// not done yet###
				[radMap initRadialMapWithTab:gaSort];
			} else {
				[traj initRadialTrajWithTab:gaTab startPoint:echoStart];
				[radMap initRadialMapWithTab:gaTab];
			}
		} else {
			[traj initRadialToshiba];			// rad180
		}

		if (thr2RR > 0) {
			[traj remove2RR:pw_c thres:thr2RR];
		}
		[traj saveAsKOImage:@"traj.img"];

[raw saveAsKOImage:@"raw1.img"];
		// rad phase correction ...
		if (goldenAngle) {
			if (sort) {
				[raw radPhaseCorr360];	// ok now ? 11-7-2017
			} else {
				[raw radPhaseCorrWithGATab:gaTab]; // ref based phase correction
			}
		} else {
			[raw radPhaseCorr];	// 0-180 counter-view
		}
TIMER_END("rad pcorr");
[raw saveAsKOImage:@"raw2.img"];
//exit(0);

// slice pre-proc
		pw = [raw copy];
		//- (RecImage *)removeSliceAtIndex:(int)ix forLoop:(RecLoop *)lp;
		[pw swapLoop:sl withLoop:pe];			// [ch, pe, *sl, *rd]
		[pw fft1d:rd direction:REC_FORWARD];	// [ch, pe, *sl,  rd]
[pw magnitude];
[pw makeComplex];

		tmp_img = [pw copy];
		[tmp_img saveAsKOImage:@"pw.img"];
		tmp_img = [tmp_img combineForLoop:ch];
	//	[tmp_img gauss2DHP:0.2];
		[tmp_img saveAsKOImage:@"pw_c.img"];

		raw = [pw copy];
		[raw swapLoop:sl withLoop:pe];
		[raw fft1d:rd direction:REC_INVERSE]; // pw -> raw [ch, sl, pe, *rd]
		[raw saveAsKOImage:@"raw3.img"];	// kx, proj, z
		[raw saveToFile:@"RAWsav" relativePath:YES];	// cache this to save time

// (3) === initial gridding
        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
//	[grid setMask:NO];	// default:YES
//	[grid setDumpRaw:YES];
		[grid setNop:nop];
        [grid grid2d:raw to:img];
		[img saveAsKOImage:@"img.img"];

        img_c = [img combineForLoop:ch];
		[img_c saveAsKOImage:@"img_c.img"];
exit(0);

//		img_coro = [img_c toCoronal];
		img_coro = [img_c copy];
		[img_coro swapLoop:[img_coro yLoop] withLoop:[img_coro zLoop]];
		[img_coro saveAsKOImage:@"img_coro.img"];
TIMER_END("initial gridding");

		img_lap = [img_coro copy];
		[img_lap laplace3d:REC_FORWARD];
		[img_lap magnitude];
	//	[img_lap gradMag3d];

		[img_lap gauss3DLP:0.1];
		[img_lap saveAsKOImage:@"img_coro_l.img"];

	// inner / outer ## merge center pos with def vector calc
		if (1) {
			// param for central mask
			float			rx = 0.32, ry = 0.22, d = 0.3, cx = 0.07, cy = -0.1;
			printf("central mask\n");
			imgo = [img copy];
			[imgo fermiWithRx:rx ry:ry d:d x:cx y:cy invert:YES half:NO];	// outer
			[imgo saveAsKOImage:@"imgo.img"];
			pwo = [imgo reprojectWithMap:radMap];
			[pwo copyLoopsOf:pw];
			[pwo saveAsKOImage:@"pwo.img"];
			pwi = [pw copy];
			[pwi subImage:pwo];
			[pwi saveAsKOImage:@"pwi.img"];
			imgi = [img copy];
			[imgi fermiWithRx:rx ry:ry d:d x:cx y:cy invert:NO half:YES];	// inner (initial recon)
		}

// (4) === iterative rigid-body correction
	//	if (rigid) {
		if (0) {
			printf("rigid body correction\n");
			sft = rigidCorr(pwi, imgi, 1);
			[sft saveAsKOImage:@"sft_out.img"];
			exit(0);
		}

// == testing 3D deformation ### === 5/31/2017, restarted on 6/27/2017
//		(rigidCorr has better result... should be the same now ######)
		if (1) {
			RecImage	*pwBins, *bTab;
			RecImage	*pwb, *rawb, *trajb;	// probably no need for pwb
			RecImage	*imgb[3], *def;
			RecImage	*sftb;
			int			i, j;
			float		mn;
			int			nBin		= 3;

			// detect 1D shift
			pwr = [imgi reprojectWithMap:radMap];					// reprojection first...
			corr = [pwi xyCorrelationWith:pwr width:0.2 triFilt:YES];	// 0.2
			sft = [corr corrToSftZ];	// pixels

			// 1D non-rigid correction ####
			[sft dumpShift];
			if (1) {
				float			scl;
				int				nscl = 5;	// 20
				RecImage		*fcs, *imgs, *mxSft;
				RecLoop			*scLp;
				RecLoopControl	*lc;

				scLp = [RecLoop loopWithDataLength:nscl];
				// coronal
				fcs  = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];
				imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];

				grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
				for (i = 0; i < nscl; i++) {
					scl = 10.0 / nscl * (i + 1);
					printf("scale %d = %f\n", i, scl);
					sftb = [sft copy];
					[sftb multByConst:scl];
					pws = [pw correctZShift:sftb];

					raw = [pws copy];
					[raw fft1d:[raw xLoop] direction:REC_INVERSE];
					[raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];

					img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
					[grid grid2d:raw to:img];
					img_c = [img combineForLoop:ch withCoil:coil];
					path = [NSString stringWithFormat:@"imgs%d.img", i];
					[img_c saveAsKOImage:path];
					[imgs copySlice:img_c atIndex:i forLoop:scLp];
					img_coro = [img_c copy];

					[img_coro laplace3d:REC_FORWARD];
					[img_coro magnitude];
					[img_coro gauss3DLP:0.1];
					[img_coro divImage:img_lap];
					[fcs copySlice:img_coro atIndex:i forLoop:scLp];
				}
				[fcs saveAsKOImage:@"imgLapScl.img"];
				[imgs saveAsKOImage:@"imgs.img"];

			// find best focus index
				mxSft = [fcs maxIndexForLoop:scLp];
				[mxSft saveAsKOImage:@"imgMxSft.img"];

			// select best image
				img_coro = [imgs selectSft:mxSft];
				[img_coro saveAsKOImage:@"img_final_coro.img"];
				[img_c copyImage:img_coro];
				[img_c saveAsKOImage:@"img_final.img"];
			}
	exit(0);

			// reconstruct each bin
			// pwBins [proj]
			pwBins = [sft binTabWithNBins:nBin binSize:[sft xDim] / 2];
			[pwBins saveAsKOImage:@"bins.img"];
			for (i = 0; i < nBin; i++) {
				bTab = [pwBins sliceAtIndex:i forLoop:[pwBins yLoop]];
				pwb = [pwi replaceLoop:pe withTab:bTab];	// 0: exp, 1: mid, 2 inh
				if (0) {	// not necessary (no improvements)
					sftb = [sft replaceLoop:pe withTab:bTab];
					path = [NSString stringWithFormat:@"sft_b%d.img", i];
					[sftb saveAsKOImage:path];
					mn = [sftb meanVal];
					for (j = 0; j < [sftb xDim]; j++) {
						[sftb data][j] -= mn;
					}
					pwb = [pwb correctZShift:sftb];
				}

				rawb = [pwb copy];
				[rawb fft1d:[rawb xLoop] direction:REC_INVERSE];
				[rawb swapLoop:sl withLoop:[pwb zLoop]];

				trajb = [traj replaceLoop:pe withTab:bTab];
				grid = [RecGridder gridderWithTrajectory:trajb andRecDim:[img xDim]];
				[grid setNop:nop];
				[grid grid2d:rawb to:img];
				path = [NSString stringWithFormat:@"imgb%d.img", i];
				[img saveAsKOImage:path];
				imgb[i] = [img combineForLoop:ch];
				path = [NSString stringWithFormat:@"imgb%d_c.img", i];
				[imgb[i] saveAsKOImage:path];
				img_coro = [RecImage imageWithImage:imgb[i]];
				[img_coro swapLoop:[imgb[i] zLoop] withLoop:[imgb[i] yLoop]];
				[img_coro copyImage:imgb[i]];
				path = [NSString stringWithFormat:@"imgb%d_coro.img", i];
				[img_coro saveAsKOImage:path];
			}
exit(0);
			// === calc 3D def vector field ===
			def = [imgb[1] defVectorWithRef:imgb[0]];
			[def saveAsKOImage:@"def.img"];	// 3 x 3 x 3 blocks

			sftb = [def defToSftWithTab:gaSort sft:sft]; // sft [blk, pe]
			[sftb saveAsKOImage:@"sft_b.img"];

			grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
			[grid setNop:nop];
			// loop over all blocks
			for (i = 0; i < 27; i++) {
				@autoreleasepool {
					printf("block %d\n", i);
					sft = [sftb sliceAtIndex:i forLoop:[sftb topLoop]];	// pixels
					pws = [pw correctShift:sft];
					tmp_img = [pws combineForLoop:ch];
					[tmp_img gauss2DHP:0.2];
					path = [NSString stringWithFormat:@"pws%d_c.img", i];
					[tmp_img saveAsKOImage:path];
					raw = [pws copy];
					[raw swapLoop:sl withLoop:pe];
					[raw fft1d:rd direction:REC_INVERSE]; // pw -> raw [ch, sl, pe, *rd]
					[grid grid2d:raw to:img];
				//	[img saveAsKOImage:@"img_s.img"];
					tmp_img = [img combineForLoop:ch];
					path = [NSString stringWithFormat:@"imgbs%d_c.img", i];
					[tmp_img saveAsKOImage:path];
				}
			}


			exit(0);
		}

	}	// autorelease pool

	return 0;
}

void
test_diff(RecImage *pw, RecImage *img)
{
	RecImage	*pwr, *err, *raw, *img_s, *img_corr;

	[img saveAsKOImage:@"img.img"];
	[pw saveAsKOImage:@"pw.img"];
	img_s = [RecImage imageWithImage:img];
	pwr = [img reprojectWithMap:radMap];
	[pwr saveAsKOImage:@"pwr.img"];
	err = [pw copy];
	[err subImage:pwr];
	[err saveAsKOImage:@"err_pw.img"];
	raw = [err pwToSin];
	[raw saveAsKOImage:@"err_sin.img"];
	[raw fft1d:rd direction:REC_INVERSE];   // sin -> raw
	// img <- raw
	[grid grid2d:raw to:img_s];	// freq/img not correct
	[img_s saveAsKOImage:@"err_img.img"];
	img_corr = [img copy];
	[img_corr addImage:img_s];
	[img_corr saveAsKOImage:@"err_corr1.img"];
	
	img_s = [img_corr combineForLoop:ch];
	[img_s saveAsKOImage:@"err_corr_c.img"];

exit(0);
	// 2nd
	pwr = [img reprojectWithMap:radMap];
	[pwr saveAsKOImage:@"pwr2.img"];
	err = [pw copy];
	[err subImage:pwr];
	[err saveAsKOImage:@"err_pw2.img"];
	raw = [err pwToSin];
	[raw saveAsKOImage:@"err_sin2.img"];
	[raw fft1d:rd direction:REC_INVERSE];   // sin -> raw
	// img <- raw
	[grid grid2d:raw to:img_s];	// freq/img not correct
	[img_s saveAsKOImage:@"err_img2.img"];
	img_corr = [img copy];
	[img_corr subImage:img_s];
	[img_corr saveAsKOImage:@"err_corr2.img"];


	exit(0);
}

// output is sft
// single step (iterative version doesn't work well)
RecImage *
rigidCorr(RecImage *pw, RecImage *img, int maxIter)
{
	RecImage	*imgb, *pwb, *rawb, *trajb, *sftb;
	RecImage	*sft;
	RecImage	*pwBins, *bTab;
	RecImage	*pwr;
	RecImage	*img_tmp, *img_coro;
	RecImage	*corr;
	int			i = 0;
	NSString	*path;
	float		mn;

	sft = nil;
	imgb = [img copy];

//	=== first sft estimation
	pwr = [imgb reprojectWithMap:radMap];
	path = [NSString stringWithFormat:@"pwr%d.img", i];
	img_tmp = [pwr combineForLoop:ch];
	[img_tmp gauss2DHP:0.2];
	[img_tmp saveAsKOImage:path];

	corr = [pw xyCorrelationWith:pwr width:0.2 triFilt:YES];	// 0.2
	path = [NSString stringWithFormat:@"corr%d.img", i];
	[corr saveAsKOImage:path];
	sft = [corr corrToSftZ];	// uses PCA, z only, unit:pixels
//	sft = [corr corrToSft];	// uses PCA, xy shift, unit:pixels
	path = [NSString stringWithFormat:@"sft%d.img", i];
	[sft saveAsKOImage:path];

// === divide into bins
	pwBins = [sft binTabWithNBins:3 binSize:[sft xDim] / 2];

	for (i = 0; i < 3; i++) {
//int sc;
//for (sc = 0; sc < 2; sc++) {
		bTab = [pwBins sliceAtIndex:i forLoop:[pwBins yLoop]];
		pwb = [pw replaceLoop:pe withTab:bTab];	// 0: exp, 1: mid, 2 inh
		img_tmp = [pwb combineForLoop:ch];
		[img_tmp gauss2DHP:0.2];
		path = [NSString stringWithFormat:@"pw_b%d.img", i];
		[img_tmp saveAsKOImage:path];
// === shift correction toward bin center
		// sft : real only, single channel [pe:256]
		sftb = [sft replaceLoop:pe withTab:bTab];
		path = [NSString stringWithFormat:@"sft_b%d.img", i];
		[sftb saveAsKOImage:path];
//		[sftb subtractMeanForLoop:[sftb xLoop]];

		if (1) {
			int j;
			mn = [sftb meanVal];
			for (j = 0; j < [sftb xDim]; j++) {
			//	printf("%d %f\n", j, [sftb data][j]);
				[sftb data][j] -= mn;
			}
		}
	// correction ####
	//	pwb = [pwb correctShift:sftb];
//[sftb multByConst:(float)sc];
		pwb = [pwb correctZShift:sftb];
		img_tmp = [pwb combineForLoop:ch];
		[img_tmp gauss2DHP:0.2];
		path = [NSString stringWithFormat:@"pws_b%d.img", i];
		[img_tmp saveAsKOImage:path];

// === recon each bin
		trajb = [traj replaceLoop:pe withTab:bTab]; //-> make traj directly
		rawb = [pwb copy];
		[rawb fft1d:[rawb xLoop] direction:REC_INVERSE];
		[rawb swapLoop:sl withLoop:[pwb zLoop]];
	//	[rawb saveAsKOImage:@"rawb.img"];
	//	[trajb saveAsKOImage:@"trajb.img"];
		grid = [RecGridder gridderWithTrajectory:trajb andRecDim:[imgb xDim]];
		[grid grid2d:rawb to:imgb];
		img_tmp = [imgb combineForLoop:ch];
		path = [NSString stringWithFormat:@"img_b%d.img", i];
		[img_tmp saveAsKOImage:path];

		img_coro = [RecImage imageWithImage:img_tmp];
		[img_coro swapLoop:[img_coro zLoop] withLoop:[img_coro yLoop]];
		[img_coro copyImage:img_tmp];
		path = [NSString stringWithFormat:@"imgb%d_coro.img", i];
		[img_coro saveAsKOImage:path];
//}
	}


// #######
exit(0);

	path = [NSString stringWithFormat:@"pwb%d.img", i];
	img_tmp = [pwb combineForLoop:ch];
	[img_tmp gauss2DHP:0.2];
	[img_tmp saveAsKOImage:path];

	trajb = [traj replaceLoop:pe withTab:bTab]; //-> make traj directly
	rawb = [pwb copy];
	[rawb fft1d:[rawb xLoop] direction:REC_INVERSE];
	[rawb swapLoop:sl withLoop:[pwb zLoop]];
	[rawb saveAsKOImage:@"rawb.img"];
	[trajb saveAsKOImage:@"trajb.img"];
	grid = [RecGridder gridderWithTrajectory:trajb andRecDim:[imgb xDim]];
	[grid grid2d:rawb to:imgb];
	img_tmp = [imgb combineForLoop:ch];
	path = [NSString stringWithFormat:@"img_s%d.img", i];
	[img_tmp saveAsKOImage:path];


	return sft;
}

