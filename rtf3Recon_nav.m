//
//	rtf3 recon (Toshiba version, + nav)
//
//	K. Oshio	12-12-2017
//
//  merge with rtf3Recon when finished
//
// === plans ===
//	*forked off from rtf3Recon_t		12-12-2017
//	*try sharper binned recon using nav data (correct each bin)	3-7-2019
//	*create reduced-view images for testing
//	parallel imaging for radial sampling -> suspend ... phantom study ?
//  merge shift est part from corTest.m - test_26()
//	test liver data (2019-12-23)
//		make switch for nav / central est
//		include central est code (just call lib)
//		
//		

#import <RecKit/RecKit.h>
#import <sys/time.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

// global
RecLoop		*ch, *sl, *pe, *rd;
RecLoop		*kx, *ky;
RecImage	*radMap, *traj;
int			coil = Coil_None;	// TOSHIBA_15;
RecGridder	*grid;

int
main(int ac, char *av[])
{
	NSString		*path;
	RecVector		fov;
	RecImage		*raw, *img, *pw, *pws, *img_c, *img_hi;
	RecImage		*gaTab = nil;		// theta tab, in time order
	RecImage		*gaSort;			// theta tab, k order, with orig index
	RecImage		*sft = nil;
	RecImage		*tmp_img;

    RecImage        *img_coro;

	int				i, nCh, nSlc;
	int				recDim = 256;
	int				zDim = 128;		// slice zero-fill

	// recon params
	BOOL			halfFT		= NO;
	BOOL			goldenAngle	= YES;	// ?? not used anymore ???
	BOOL			navData		= NO;
	int				echoStart	= 0; //112;
	NSString		*base = @"../toshiba_images";
	NSString		*date, *run;
	int				ser;
	RecLoop			*subPe;

TIMER_ST
    @autoreleasepool {
        system("rm *.img");
        system("rm IMG_*");
        system("rm sft*.txt");

	//	[RecImage setFFTdbg:NO];	// turn off FFT unit check

// (1) === load raw data
		printf("1: load raw data\n");
	//	date = @"2017-10-24"; run = @"/Run23660.-5213.15";
	//	date = @"2018-02-19"; run = @"Run266.-5003.15";			// heart + nav
		date = @"2019-12-23"; run = @"Run65418.-5004.04-PW";	// liver (no nav)
		ser = 15;
		path = [NSString stringWithFormat:@"%@/%@/%@", base, date, run];
		raw = [RecImage imageWithToshibaFile:path vorder:&gaTab fov:&fov];
		if (navData) {
			path = [NSString stringWithFormat:@"%@/%@/probe%d.txt", base, date, ser];
			sft = Rec_load_edge(path);
		}
		echoStart = 0;	// SSFP
		halfFT = NO;

		if (raw == nil) {
			printf("file not found [%s]\n", [path UTF8String]);
			exit(0);
		}
	// raw data loops
		rd = [RecLoop findLoop:@"Read"];
		pe = [RecLoop findLoop:@"Phase"];
		ch = [RecLoop findLoop:@"Channel"];
		sl = [RecLoop findLoop:@"Slice"];
		kx = [RecLoop loopWithName:@"kx" dataLength:recDim];
		ky = [RecLoop loopWithName:@"ky" dataLength:recDim];
		nCh = [ch dataLength];
		nSlc = [sl dataLength];

	// crop (reduce number of projections for parallel imaging test)
	// ###### gaTab, raw, sft
	if (0) {
		subPe = [RecLoop loopWithDataLength:[pe dataLength]/2];	// 2
		[gaTab replaceLoop:[gaTab xLoop] withLoop:subPe];
		[raw replaceLoop:pe withLoop:subPe];
		[sft replaceLoop:[sft xLoop] withLoop:subPe];
		pe = subPe;
	}
	printf("pe dim = %d\n", [pe dataLength]);

		gaSort = Rec_goldenAngleSort(gaTab);	// gaTab is theta tab passed by xml header
		[raw saveAsKOImage:@"raw0.img"];		// gaSort is sorted version
		[gaSort saveAsKOImage:@"IMG_gasort"];
		[gaTab saveAsKOImage:@"IMG_gatab"];

// ### compare
		if (0) {
			tmp_img = [RecImage imageWithKOImage:@"single_sav/IMG_sft1.img"];
			[sft copyImageData:tmp_img];
			[sft negate];
		}

TIMER_END("load raw");

		if (halfFT) {
			sl = [raw halfFT2:sl];				// [ch, pe,  sl,  rd] (pw)
		} else {
			sl = [raw sftCenter:sl zeroFillTo:zDim];
		}
	
// (2) === pre-processing
		printf("2: pre-processing\n");
		// create radial trajectory
		traj  = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rd, nil];
		radMap = [RecImage imageOfType:RECIMAGE_MAP  withLoops:pe, rd, nil];

		// reorder both raw / sft using gaSort tab
		[raw reorder:pe withTab:gaSort];			// view angle order

		[traj initRadialTrajWithSortTab:gaSort startPoint:echoStart];
		[radMap initRadialMapWithTab:gaSort];
		[traj saveAsKOImage:@"traj.img"];

		// rad phase correction ...
		if (goldenAngle) {
			[raw radPhaseCorr360];	// 0-th (for error due to x-shift)
		} else {
			[raw radPhaseCorr];	// 0-180 counter-view
		}
		[raw saveAsKOImage:@"raw2.img"];	// after phase correction
	//	[raw saveToFile:@"raw2" relativePath:YES];	// cache this to save time

// slice pre-proc
		pw = [raw copy];
		[pw swapLoop:sl withLoop:pe];			// [ch, pe, *sl, *rd]
		[pw fft1d:rd direction:REC_FORWARD];	// [ch, pe, *sl,  rd]

	// 0-th order corr added...
	// higher order correction might be necessary ...
		if (0) {
			for (i = 0; i < [ch dataLength]; i++) {
				RecImage	*chImg = [pw sliceAtIndex:i forLoop:ch];
				[chImg pcorr0EachSlice];
				[pw copySlice:chImg atIndex:i forLoop:ch];
			}
		}
		[pw saveAsKOImage:@"pw.img"];

// after 0th order (image domain) phase correction
		raw = [pw copy];
		[raw swapLoop:sl withLoop:pe];
		[raw fft1d:rd direction:REC_INVERSE]; // pw -> raw [ch, sl, pe, *rd]
		[raw saveAsKOImage:@"raw3.img"];	// kx, proj, z

// (3) === initial gridding
		printf("3: initial gridding\n");
        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
		// nop = 2 is better, but much slower
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:YES];
		printf("pix size: x:%f, y:%f, z:%f\n", fov.x / [kx dataLength], fov.y / [ky dataLength], fov.z / [sl dataLength]);
		printf("recDim = %d\n", [img xDim]);

		// dbg code...
		[grid setMask:NO];
		[grid setDumpRaw:NO];

        [grid grid2d:raw to:img];
		[img saveAsKOImage:@"img.img"];
		
		[img saveToFile:@"img_sav.recimg" relativePath:YES];

        img_c = [img combineForLoop:ch];
		[img_c saveAsKOImage:@"img_c.img"];
		img_hi = [img_c oversample];
		[img_hi saveAsKOImage:@"img_c_hi.img"];
		img_coro = [img_c copy];
		[img_coro swapLoop:[img_coro yLoop] withLoop:[img_coro zLoop]];
		[img_coro saveAsKOImage:@"img_coro.img"];
TIMER_END("initial gridding");

// ### nav data
		if (navData) {
			[sft replaceLoop:[sft xLoop] withLoop:pe];
			[sft saveAsKOImage:@"nav_sft_t.img"];		// time order, 0-mean, unit is mm (up is +)
			[sft reorder:[sft xLoop] withTab:gaSort];	// view angle order
			[sft saveAsKOImage:@"nav_sft.img"];
		} else {
			// central view est ###
			// input is pw
			sft = [pw shiftFromK0];	// shift estimation from center-of-kspace view (POCS)
		}

// === test single shift -> OK
		if (0) {
			// ## chk scale
			float *p = [sft data];
			for (i = 0; i < 10; i++) {
				printf("%d %f\n", i, p[i]);
			}

			pws = [pw correctZShift:sft];
			grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:YES];
			raw = [pws copy];
			[raw fft1d:[raw xLoop] direction:REC_INVERSE];
			[raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];
			[grid grid2d:raw to:img];
			img_c = [img combineForLoop:ch];
			[img_c saveAsKOImage:@"IMG_s"];
			img_coro = [img_c copy];
			[img_coro swapLoop:[img_coro yLoop] withLoop:[img_coro zLoop]];
			[img_coro saveAsKOImage:@"imgs_coro.img"];

			exit(0);
		}

// (4) === step-shift correction (z only) ====
		if (1) {
			float			scl;
			int				nscl = 10;
			float			scl_range = 3.0;	// (0 - 2.0)
			RecImage		*imgs, *mxSft;
			RecImage		*sft_scl, *pws;
			RecLoop			*scLp;
			float			lap_w = 0.2;	// 0.3

			printf("4: step-shift\n");
			grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim] nop:1 densityCorrection:YES];

			scLp = [RecLoop loopWithDataLength:nscl];
			imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:scLp, [img yLoop], [img zLoop], [img xLoop], nil];	// coro
		//	imgs = [RecImage imageOfType:RECIMAGE_REAL withLoops:[img yLoop], scLp, [img zLoop], [img xLoop], nil];	// coro, scale first (doesn't work)

			for (i = 0; i < nscl; i++) {
				scl = (float)i * scl_range / nscl;
				printf("scale %d = %f\n", i, scl);

				sft_scl = [sft copy];			// pixels
				[sft_scl multByConst:scl];
				pws = [pw correctZShift:sft_scl]; // z shift

				tmp_img = [pws copy];
				raw = [pws copy];
				[raw fft1d:[raw xLoop] direction:REC_INVERSE];
				[raw swapLoop:[raw yLoop] withLoop:[raw zLoop]];
				[grid grid2d:raw to:img];
				img_c = [img combineForLoop:ch];
				[imgs copySlice:img_c atIndex:i forLoop:scLp];
			}
			tmp_img = [imgs copy];
			[tmp_img swapLoop:scLp withLoop:[img yLoop]];
			[tmp_img saveAsKOImage:@"imgs_scl.sav"];

// (5) === find best focus scale ====
			printf("5: select best shift\n");
			tmp_img = [imgs copy];
//			[tmp_img laplace2d:REC_FORWARD];
			[tmp_img grad1dForLoop:[tmp_img yLoop]];
			[tmp_img magnitude];
			[tmp_img square];
			[tmp_img gauss3DLP:lap_w];
			[tmp_img swapLoop:scLp withLoop:[img yLoop]];	// ## chk
			[tmp_img saveAsKOImage:@"imgs_scl_lap.img"];

			mxSft = [tmp_img peakIndexForLoop:scLp];	// sharper with larger laplacian
			[mxSft saveAsKOImage:@"imgMxSft.img"];		// ok
			[mxSft gauss2DLP:0.1];
			[mxSft saveAsKOImage:@"imgMxSft_f.img"];
			tmp_img = [imgs selectSftF:mxSft];
			[tmp_img saveAsKOImage:@"imgs_final.img"];
			[tmp_img swapLoop:[tmp_img yLoop] withLoop:[tmp_img zLoop]];
			[tmp_img saveAsKOImage:@"imgs_final_ax.img"];
		} // === focused rigid-body correction === 
	}	// autorelease pool

	return 0;
}

