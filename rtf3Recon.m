//
//	rtf3 recon (main, GE)
//
//	K. Oshio	1-27-2015
//
//  FOV / slthk / nslice is assumed to be 300, 2, 64
//  original pixel size is: xy:1.172, z:2, z/x ratio:1.707
//
//  forked off from rtf3Recon4.m    4-28-2020

//  started update for central line nav 5-2-2020
//  cmd-line arg -> hard-coded path + zflip etc


/* 2-25-2015

P16384.7	sakai
P10240.7	akasaka

LSB first
RDBM header version = 11.00
 time: 12.735471 (sec) img
====== reproject / correction 0 ======
 time: 19.664232 (sec) pwr
0 1.000000
 time: 11.822875 (sec) img_s
====== reproject / correction 1 ======
 time: 0.723146 (sec) pwr
1 0.219740
 time: 11.768156 (sec) img_s
====== reproject / correction 2 ======
 time: 0.724449 (sec) pwr
2 0.078888
====== scaling ======
 time: 52.026984 (sec) scale/shift
====== final gridding ======
 time: 100.056864 (sec) grid
 time: 218.458848 (sec) TOTAL

*/

#import <RecKit/RecKit.h>
#import <sys/time.h>
#import "RecImagePW.h"
#import <RecKit/timer_macros.h>

// volunteers
//NSString *base = @"/Users/oshio/epic/rtf3_data/volunteer";
//NSString *name = @"belt-off"; int pNum = 34304; BOOL zFlip = YES;  // signal X
//NSString *name = @"belt-on";  int pNum = 33792; BOOL zFlip = YES;  // signal X
//NSString *name = @"hashimoto1"; int pNum = 35328; BOOL zFlip = YES;  // 0.6
//NSString *name = @"hashimoto2"; int pNum = 36352; BOOL zFlip = YES;  // 0.0, 2
//NSString *name = @"okabe2-1"; int pNum = 23040; BOOL zFlip = NO;    // 0.1, 0
//NSString *name = @"okabe2-2"; int pNum = 23552; BOOL zFlip = NO;  // 0.5, 0
//NSString *name = @"okabe2-3"; int pNum = 24064; BOOL zFlip = NO;  // 0.3, 9
//NSString *name = @"okuda-nitrol"; int pNum = 14336; BOOL zFlip = NO;  // 0.2, 8

// patients
NSString *base = @"/Users/oshio/epic/rtf3_data/clinical";
//NSString *name = @"akasaka";      int pNum = 10240; BOOL zFlip = NO;  // 0.3, 5
//NSString *name = @"akita";        int pNum = 25088; BOOL zFlip = NO;  // 0.6, 7 
//NSString *name = @"akiyama";      int pNum = 24064; BOOL zFlip = NO;  // 0.7, 0
//NSString *name = @"anai";         int pNum = 13312; BOOL zFlip = YES; // 0.8, 7
//NSString *name = @"arimoto";      int pNum =  8704; BOOL zFlip = NO;  // 0.8, 3
//NSString *name = @"eno";          int pNum = 12288; BOOL zFlip = NO;  // 0.8, 3
//NSString *name = @"funatogawa";   int pNum =  8704; BOOL zFlip = NO;  // 0.6, 4
//NSString *name = @"furuyama";     int pNum = 41984; BOOL zFlip = NO;  // 0.9, 6 ox
//NSString *name = @"hatakeyama";   int pNum = 49152; BOOL zFlip = YES; // 0.4, 6
//NSString *name = @"hayashi";      int pNum = 40960; BOOL zFlip = NO;  // 0.6, 0
//NSString *name = @"ikeda";        int pNum = 13824; BOOL zFlip = YES; // 0.3. 2 -> minErr at 0.6
//NSString *name = @"ishii";        int pNum = 19456; BOOL zFlip = NO;  // 0.8, 3
//NSString *name = @"itou";         int pNum = 34304; BOOL zFlip = YES; // 0.9, 1
//NSString *name = @"kanai";        int pNum = 41984; BOOL zFlip = NO;  // 0.8, 8
//NSString *name = @"kariya";       int pNum =  3584; BOOL zFlip = NO;  // 0.9, 1
//NSString *name = @"kawano";       int pNum = 18944; BOOL zFlip = YES; // 0.8, 2
//NSString *name = @"kimura";       int pNum = 28672; BOOL zFlip = YES; // 0.8, 7
//NSString *name = @"kiura";        int pNum = 14336; BOOL zFlip = NO;      // 0.9, 6 -> X
NSString *name = @"kobayashi";    int pNum = 17920; BOOL zFlip = NO;      // 0.8, 9
//NSString *name = @"";             int pNum = ; BOOL zFlip = NO;
//NSString *name = @"";             int pNum = ; BOOL zFlip = NO;
//NSString *name = @"";             int pNum = ; BOOL zFlip = NO;
//NSString *name = @"";             int pNum = ; BOOL zFlip = NO;
//NSString *name = @"";             int pNum = ; BOOL zFlip = NO;

int
//main(int ac, char *av[])
main()
{
	NSString		*path;
	RECCTL			ctl;
	RecImage		*raw, *img, *pw;
    RecImage        *sft;
    RecImage        *img_c, *pw_c, *pwr, *pws, *pws_c;
    RecImage        *img_hi, *img_coro, *img_coro_c;
	RecLoop			*ch, *sl, *pe, *rd;
	RecLoop			*rdZF;
	RecLoop			*kx, *ky;
	RecImage		*radMap, *traj, *thTab;
    RecImage        *imgs;
    RecGridder      *grid;
    int             i;
	int				coil = GE_Card_8_new; // try image based est
 
TIMER_ST
    @autoreleasepool {
        
    system("rm *.img, img*, IMG_*");
    system("rm sft*.txt");

        if (0) {
            path = [NSString stringWithFormat:@"%@/%@/pw_sav", base, name];
            pw = [RecImage imageFromFile:path relativePath:NO];
            sft = [pw shiftFromK0];
            [sft saveAsKOImage:@"IMG_sft"];
            exit(0);
        }


	// (1) === load raw data
        path = [NSString stringWithFormat:@"%@/%@/P%05d.7", base, name, pNum];
        raw = [RecImage imageWithPfile:path RECCTL:&ctl];	// 3D
        if (raw == nil) {
            printf("Raw data not found\n");
            exit(-1);
        }
	// raw data loops
        rd = [RecLoop findLoop:@"Read"];
        pe = [RecLoop findLoop:@"Phase"];
        ch = [RecLoop findLoop:@"Channel"];
        sl = [RecLoop findLoop:@"Slice"];
		[raw removePointLoops];

    // rec dim
        rdZF = [RecLoop loopWithName:@"rd_zf" dataLength:256];
        kx   = [RecLoop loopWithName:@"kx"    dataLength:ctl.rc_xres];
        ky   = [RecLoop loopWithName:@"ky"    dataLength:ctl.rc_yres];

	// (2) === pre-processing
    // z-unzip -> handle half fourier (just shift to center, 10240 etc)
        [raw fft1d:sl direction:REC_INVERSE];
        [raw shift1d:sl];
        [raw xFlip];    // if GE... seq was wrong...
        if (zFlip) {
            [raw flipForLoop:sl];
        }

    // read-axis zero-fill (224 -> 256)
        [raw replaceLoop:rd withLoop:rdZF];
        [raw radPhaseCorr]; // radial phase correction
		[raw saveAsKOImage:@"raw.img"];

        pw = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, pe, sl, rdZF, nil];

		[raw fft1d:sl direction:REC_FORWARD];
		[raw fft1d:rdZF direction:REC_FORWARD];	// sinogram
        [pw copyImage:raw]; // sinogram -> pw
		pw_c = [pw combinePWForLoop:ch withCoil:coil];
        [pw_c saveAsKOImage:@"pw_c.img"];
        [pw saveAsKOImage:@"pw.img"];

        // save to original data location
        path = [NSString stringWithFormat:@"%@/%@/pw_sav", base, name];
        [pw saveToFile:path relativePath:NO];
   
        traj = [RecImage imageOfType:RECIMAGE_KTRAJ withLoops:pe, rdZF, nil];

        thTab = [traj initRadialTraj];          // for gridding
        radMap = [traj mapForRadial];   // for reprojection
        [raw fft1d:rdZF direction:REC_INVERSE]; // pw -> raw

	// (3) === initial gridding
        img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ch, sl, ky, kx, nil];
        grid = [RecGridder gridderWithTrajectory:traj andRecDim:[img xDim]];
        [grid grid2d:raw to:img];
        img_c = [img combineForLoop:ch withCoil:coil];
		[img_c SCIC];
		[img_c saveAsKOImage:@"img_c.img"]; // first image (before correction)

		img_coro_c = [RecImage imageWithImage:img_c];
		[img_coro_c swapLoop:[img_c zLoop] withLoop:[img_c yLoop]];
		[img_coro_c copyImage:img_c];
		[img_coro_c saveAsKOImage:@"img_coro.img"];

    // sft est
        sft = [pw shiftFromK0];
        [sft saveAsKOImage:@"IMG_sft"];

    // correction
        imgs = [img stepCorrWithPW:pw gridder:grid sft:sft];
        [imgs saveAsKOImage:@"IMG_imgs"];


    } // autoreleasepool
TIMER_TOTAL

	return 0;
}
