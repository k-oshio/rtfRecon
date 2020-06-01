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

//  plans:
//  90% of cases are ok with global initial condicion search (ref slice, st frac)
//  *remaining cases have large deformation -> correction won't work anyway
//  chk failing case with apparently benign resp pattern (furuyama, sakai, shinomiya)
//  *oscillation in y ? try damping
//  wider window


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
//NSString *name = @"akasaka";      int pNum = 10240; BOOL zFlip = NO;  // 0.7, 7 X (bad data)
//NSString *name = @"akita";        int pNum = 25088; BOOL zFlip = NO;  // 0.9, 1 X (bad data)
//NSString *name = @"akiyama";      int pNum = 24064; BOOL zFlip = NO;  // 0.7, 0 X -> info in upper part
//NSString *name = @"anai";         int pNum = 13312; BOOL zFlip = YES; // 0.9, 3
//NSString *name = @"arimoto";      int pNum =  8704; BOOL zFlip = NO;  // 0.3, 3
//NSString *name = @"eno";          int pNum = 12288; BOOL zFlip = NO;  // 0.8, 2 X (ok with HPF, 0.2, 0.7, tri on, zf off)
//NSString *name = @"funatogawa";   int pNum =  8704; BOOL zFlip = NO;  // 0.1, 5
//NSString *name = @"furuyama";     int pNum = 41984; BOOL zFlip = NO;  // 0.5, 5 X -> chk
//NSString *name = @"hatakeyama";   int pNum = 49152; BOOL zFlip = YES; // 0.0, 8
//NSString *name = @"hayashi";      int pNum = 40960; BOOL zFlip = NO;  // 0.7, 9
//NSString *name = @"ikeda";        int pNum = 13824; BOOL zFlip = YES; // 0.0. 2
//NSString *name = @"ishii";        int pNum = 19456; BOOL zFlip = NO;  // 0.7, 5
//NSString *name = @"itou";         int pNum = 34304; BOOL zFlip = YES; // 0.0, 4
//NSString *name = @"kanai";        int pNum = 41984; BOOL zFlip = NO;  // 0.0, 2
//NSString *name = @"kariya";       int pNum =  3584; BOOL zFlip = NO;  // 0.8, 1
//NSString *name = @"kawano";       int pNum = 18944; BOOL zFlip = YES; // 0.9, 9
//NSString *name = @"kimura";       int pNum = 28672; BOOL zFlip = YES; // 0.9, 1
//NSString *name = @"kiura";        int pNum = 14336; BOOL zFlip = NO;  // 0.9, 7 X
//NSString *name = @"kobayashi";    int pNum = 17920; BOOL zFlip = NO;  // 0.8, 4
//NSString *name = @"kumahara";     int pNum = 17920; BOOL zFlip = NO;  // 0.9, 1
//NSString *name = @"masuda";       int pNum = 13824; BOOL zFlip = NO;  // 0.9, 9 o (focus on spleen)
//NSString *name = @"matsuda";      int pNum = 20992; BOOL zFlip = YES;  // 0.5, 8
//NSString *name = @"miyazawa_nb";  int pNum = 5632; BOOL zFlip = NO;    // 0.1, 9
//NSString *name = @"nakazawa";     int pNum = 27648; BOOL zFlip = NO;   // 0.3, 9
//NSString *name = @"noguchi";      int pNum = 36864; BOOL zFlip = NO;   // 0.8, 9
//NSString *name = @"ogino";        int pNum = 46592; BOOL zFlip = YES;  // 0.1, 6
//NSString *name = @"oka";          int pNum = 12800; BOOL zFlip = NO;   // 0.7, 1 X
//NSString *name = @"saitou";       int pNum = 15872; BOOL zFlip = YES;  // 0.9, 5
NSString *name = @"sakai";        int pNum = 16384; BOOL zFlip = YES;  // 0.9, 6 oo
//NSString *name = @"shimojo";      int pNum = 40960; BOOL zFlip = NO;   // 0.9, 9
//NSString *name = @"shinomiya";    int pNum = 13824; BOOL zFlip = YES;   // 0.9, 3
//NSString *name = @"suzuki";       int pNum = 20992; BOOL zFlip = YES;   // 0.4, 9
//NSString *name = @"takada";       int pNum = 12288; BOOL zFlip = YES;   // 0.1, 1
//NSString *name = @"tanaka";       int pNum = 32768; BOOL zFlip = YES;   // 0.5, 0
//NSString *name = @"taniguchi";    int pNum = 13824; BOOL zFlip = YES;   // 0.2, 4
//NSString *name = @"tsurufuji";    int pNum = 25088; BOOL zFlip = YES;   // 0.3, 5 X ok with ZF
//NSString *name = @"wada";         int pNum = 19456; BOOL zFlip = YES;   // 0.9, 6 X
//NSString *name = @"yamamoto";     int pNum = 43520; BOOL zFlip = NO;    // 0.9, 2 X ok with ZF or tri off
//NSString *name = @"yasue";        int pNum = 9216; BOOL zFlip = YES;    // 0.3, 0
//NSString *name = @"yawata";       int pNum = 1024; BOOL zFlip = NO;     // 0.1, 3
//NSString *name = @"yoshinari";    int pNum = 13824; BOOL zFlip = NO;    // 0.1, 9 X ok with neg damping (1.1)

// current best comb: HPF 0.8, cos filt 20, zerofill off

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
        
    system("rm *.img img* IMG_*");

        if (1) {
            path = [NSString stringWithFormat:@"%@/%@/pw_sav", base, name];
            pw = [RecImage imageFromFile:path relativePath:NO];
            sft = [pw shiftFromK0];
            [sft saveAsKOImage:@"IMG_sft"];
            pw = [pw correctZShift:sft];
            ch = [RecLoop findLoop:@"Channel"];
            pw = [pw combineForLoop:ch];
            [pw saveAsKOImage:@"IMG_pws"];
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
        [pw_c saveAsKOImage:@"pw.img"];
//        [pw saveAsKOImage:@"pw.img"];

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
        path = [NSString stringWithFormat:@"%@/%@/img.img", base, name];
		[img_c saveAsKOImage:path]; // first image (before correction)

		img_coro_c = [RecImage imageWithImage:img_c];
		[img_coro_c swapLoop:[img_c zLoop] withLoop:[img_c yLoop]];
		[img_coro_c copyImage:img_c];
        path = [NSString stringWithFormat:@"%@/%@/img_coro.img", base, name];
		[img_coro_c saveAsKOImage:path];

    // sft est
        sft = [pw shiftFromK0];
    //    [sft saveAsKOImage:@"IMG_sft"];
    
    // rigid body correction with scale = 1
        pws = [pw correctZShift:sft];
        pws = [pws combineForLoop:ch];
        [pws saveAsKOImage:@"IMG_pws"];
    
    // correction
        imgs = [img stepCorrWithPW:pw gridder:grid sft:sft];
        path = [NSString stringWithFormat:@"%@/%@/img_scl.img", base, name];
        [imgs saveAsKOImage:path];


    } // autoreleasepool
TIMER_TOTAL

	return 0;
}
