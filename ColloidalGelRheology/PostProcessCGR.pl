#! /usr/bin/perl
use Math::Trig;
$data = $ARGV[0];
$laying = 0;
open (IN_data, "< ${data}");
#$output = "${input}" . ${n}. '.pov';

############################
# Bond parameters need to be set
$kn = 20;
$ks = 20;
$kb = 20;
$kn3 = 80000;
############################

# Layer number
$y1_walls = 1;
$y2_particles = 2;
$y3_bond = 3;
$y4_bond_normal = 4;
$y5_bond_sliding = 5;
$y6_bond_bending = 6;

$i = index($data, '_W', 0)+2;
$j = index($data, 'H', $i-1);
$L = substr($data, $i, $j-$i);

$L2 = $L/2;
print "$data\n";
$checkdim = index($data, '_2D_');
printf "$checkdim\n";

if ($checkdim == -1){
    $dim = 3;
    printf "3D\n";
} else {
    $dim = 2;
    printf "2D\n";
}

$i = index($data, 'conf_', 0)+5;
$j = index($data, '.dat', $i-1);
$name = substr($data, $i, $j-$i);

printf "L = $L\n";
$output = "$name.yap";
$output2 = "stress_$name.dat";
$output3 = "eq_$name.dat";
$output4 = "energy_$name.dat";

open (OUT, "> ${output}");
open (OUT2, "> ${output2}");
open (OUT3, "> ${output3}");
open (OUT4, "> ${output4}");

# &printHead;
&yaplotColor;

$first = 1;
$initial = 1;
$first_output = 1;
$p0 = 0;
$onetime = 1;

sub spherepart {
    pi*(4.0/3 + $lackz*$lackz*$lackz/3 - $lackz*$lackz);
}
$cnt = 0;
$out2_first = 1;

while (1){
	$line = <IN_data>;
    ($buf, $equiv, $vf, $shear_strain,
	$particpe_pressure, $shear_stress,
	$lx, $lz, $rn, $rs, $rb, $rt, $ac) = split(/\s+/, $line);
	##
    if ($initial == 1){
        $equiv = "e";
        $initial = 0;
    }
	$w0z = -$lz/2;
	$w1z = $lz/2;
    if ($equiv eq "e"){
        printf  "equi\n";
        printf OUT3 "1\n";
    } else {
        printf OUT3 "0\n";
    }

	$equiv = "e";
	# import Particle position
    $line = <IN_data>;
    ($buf, $np) = split(/\s+/, $line);
    if ($buf ne "P"){last; }
	
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_data> ;
        ($x, $y, $z, $q0, $q1, $q2, $q3, $ic,$wg, $wl, $cn) = split(/\s+/, $line);    
        $wlvec[$i]=$wg;
        if ($yoko == 1){
            $posx[$i] = $z; $posy[$i] = $y; $posz[$i] = -$x;
        } else {
            $posx[$i] = $x; $posy[$i] = $y; $posz[$i] = $z;
        }
    }
    $line = <IN_data> ;
    ($buf, $nb ) = split(/\s+/, $line);
    if ($buf ne "B"){last; }
    for ($i = 0; $i < $np; $i ++){
        $stress[$i] =0;
    }

    for ($i = 0; $i < $nb; $i ++){
        $line = <IN_data> ;
        ($p0, $p1, $ib, $st, $fn, $fs, $mb, $mt) = split(/\s+/, $line);
        if ($p0 eq "#" ){
            # This is for the case that
            # nb is not correct.
            printf "bond number\n";
            last;
        }        
        ($x0[$i], $y0[$i], $z0[$i]) = ($posx[$p0], $posy[$p0], $posz[$p0]);
        ($x1[$i], $y1[$i], $z1[$i]) = ($posx[$p1], $posy[$p1], $posz[$p1]);
        $b_fn[$i] = $fn;
        $b_fs[$i] = $fs;
        $b_mb[$i] = $mb;
        $b_mt[$i] = $mt;
		$initialbond[$i] = $ib;
        $stress[$p0] += 0.5*(abs($fn) + abs($fs));
        $stress[$p1] += 0.5*(abs($fn) + abs($fs));
    }
	
	#	for ($k=0; $k < 20; $k++ ){
	if ($equiv eq "e"){
		#for ($k=0; $k < 20; $k++ ){
		if ($first_output == 1) {
			$first_output = 0;
		} else {
			printf OUT "\n" ;
		}
		&modifyCoordinateForPeriodicBoundary;
		&outputYaplot;
		$ene_n_total = 0;
		$ene_s_total = 0;
		$ene_b_total = 0;
		for ($i = 0 ; $i < $nb; $i++){
			$dist = sqrt(($x0[$i] - $x1[$i])**2 + ($y0[$i]-$y1[$i])**2 + ($z0[$i]-$z1[$i])**2);
			if ($dist > 2){
				$ene_n = 0.5*($b_fn[$i])**2 / $kn;
			} else {
				$gap = $dist - 2;
				$ene_n = 0.5*$kn*($gap)**2 + 0.25*$kn3*($gap)**4;
			}
			$ene_s = 0.5*($b_fs[$i])**2 / $ks;
			$ene_b = 0.5*($b_mb[$i])**2 / $kb;
			$ene_n_total += $ene_n;
			$ene_b_total += $ene_b;
			$ene_s_total += $ene_s;
		}
		#		if ($ene_n_total < 1e-10){
		#		$ene_n_total = 0;
		#}
		#if ($ene_s_total < 1e-10){
		#		$ene_s_total = 0;
		#}
		#if ($ene_b_total < 1e-10){
		#		$ene_b_total = 0;
		#}
		$ene_n_total = $ene_n_total / $nb;
		$ene_s_total = $ene_s_total / $nb;
		$ene_b_total = $ene_b_total / $nb;
		
		printf OUT4 "$vf $ene_n_total $ene_s_total $ene_b_total $nb \n";
		#}
	} else {
		if ($first == 1){
			$first = 0;
		} else {
			#	printf OUT "\n" ;
		}
		#&outputYaplot;
	}
		#last if ($equiv eq "n");
	#}
}
close (OUT);
close (OUT2);
close (OUT3);

sub outputYaplot
{
	&outputWalls;
	printf OUT "y $y2_particles \n";
	printf OUT "r 1\n";
	printf OUT "@ 2\n";
	
	if ($out2_first  == 1){
		printf OUT2 "$np\n";
		$out2_first = 0;
	}
	printf OUT "@ 2\n";
	for ($i = 0; $i < $np; $i ++){
		#            if ($wlvec[$i] == 1){
		#  printf OUT "@ 3\n";
		#      #Particles touching to walls
		#} elsif ($wlvec[$i] == 2){
		#   printf OUT "@ 4\n";
		#} elsif ($wlvec[$i] == -1){
		#}
		printf OUT "c $posx[$i] $posy[$i] $posz[$i]\n";
		printf OUT2 "$posx[$i] $posy[$i] $posz[$i] $stress[$i]\n";
	}
	
	printf OUT "@ 3\n";
	printf OUT "y $y4_bond_normal \n" ;
	for ($i = 0; $i < $nb; $i ++){
		$force = $b_fn[$i];
		#		printf OUT "r $force \n" ;
		if (1|| $initialbond[$i] == 1 ){
			$y0[$i] = -0.1;
			$y1[$i] = -0.1;
			printf OUT "l $x0[$i] $y0[$i] $z0[$i] $x1[$i] $y1[$i] $z1[$i]\n";
			#						if ($force > 0 ){
			#				printf OUT "@ 5 \n" ;
			#				printf OUT "r $force \n" ;
			#				printf OUT "s $x0[$i] $y0[$i] $z0[$i] $x1[$i] $y1[$i] $z1[$i]\n";
			#			} else {
			#				printf OUT "@ 3 \n" ;
			#				printf OUT "r -$force \n" ;
			#				printf OUT "s $x0[$i] $y0[$i] $z0[$i] $x1[$i] $y1[$i] $z1[$i]\n";
			#			}

		}
		
	}
}

sub modifyCoordinateForPeriodicBoundary
{
    for ($i = 0; $i < $nb; $i ++){
        if ( $laying == 1 ){
            if ( abs($z0[$i] - $z1[$i]) > 5){
                if( $z0[$i] > $z1[$i] ){
                    $z0[$i] = $z0[$i] - $L;
                } else {
                    $z1[$i] = $z1[$i] - $L;
                }
            }
            if ( abs($y0[$i] - $y1[$i]) > 5){
                if( $y0[$i] > $y1[$i] ){
                    $y0[$i] = $y0[$i] - $L;
                } else {
                    $y1[$i] = $y1[$i] - $L;
                }
            }
        } else {
            if ( abs($x0[$i] - $x1[$i]) > 5){
                if( $x0[$i] > $x1[$i] ){
                    $x0[$i] = $x0[$i] - $lx;
                } else {
                    $x1[$i] = $x1[$i] - $lx;
                }
            }
            if ( abs($y0[$i] - $y1[$i]) > 5){
                if( $y0[$i] > $y1[$i] ){
                    $y0[$i] = $y0[$i] - $lx;
                } else {
                    $y1[$i] = $y1[$i] - $lx;
                }
            }
			
			if ( abs($z0[$i] - $z1[$i]) > 5){
                if( $z0[$i] > $z1[$i] ){
                    $z0[$i] = $z0[$i] - $lz;
                } else {
                    $z1[$i] = $z1[$i] - $lz;
                }
            }
        }
    }
}

sub outputWalls
{
    if ($laying == 1){

        printf OUT "y $y1_walls \n";
        printf OUT "@ 2\n";
		#        $textposition = $w1z +  5;
        #printf OUT "t $textposition 0 0 vf=$vf \n";
        printf OUT "t 50 0 10 vf=$vf \n";
        if ($dim == 3){
            printf OUT "l $w0z -$L2 -$L2 $w0z  $L2 -$L2 \n";
            printf OUT "l $w0z  $L2 -$L2 $w0z  $L2  $L2 \n";
            printf OUT "l $w0z  $L2  $L2 $w0z -$L2  $L2 \n";
            printf OUT "l $w0z -$L2  $L2 $w0z -$L2 -$L2 \n";
            printf OUT "l $w0z  $L2  0  $w0z -$L2  0 \n";
        }
        printf OUT "l $w0z  0   $L2 $w0z 0  -$L2 \n";
        if ($dim == 3){
            printf OUT "l $w1z -$L2 -$L2 $w1z  $L2 -$L2 \n";
            printf OUT "l $w1z  $L2 -$L2 $w1z  $L2  $L2 \n";
            printf OUT "l $w1z  $L2  $L2 $w1z -$L2  $L2 \n";
            printf OUT "l $w1z -$L2  $L2 $w1z -$L2 -$L2 \n";
            printf OUT "l $w1z  $L2  0  $w1z -$L2  0 \n";
        }
        printf OUT "l $w1z  0   $L2 $w1z 0  -$L2 \n";            
        #printf OUT "y 2\n";
        #printf OUT "@ 2\n";
    } else {
        $textposition1 = $w1z +  5;
        $textposition2 = $w1z +  3;
        $textposition3 = $w1z +  1;
        printf OUT "y $y1_walls \n";
        printf OUT "@ 0\n";
        printf OUT "r 0.1\n";

		#        printf OUT "t -5 0 $textposition1 vf=$vf \n";
        #printf OUT "t -5 0 $textposition2 Strain=$strain_x  \n";
		#        printf OUT "t -5 0 $textposition3 Stress=$stress_x  \n";
		$t1=$L2+6;
		$t2=$L2+2;
		#printf OUT "@ 2\n";
		#printf OUT "t $t1 0 20 Area fraction=$vf \n";
		#printf OUT "t $t2 0 20 Compressive stress=$stress_z  \n";
		printf OUT "@ 0\n";
		#printf OUT "l -$L2 0 0  $L2 0 0 \n";
		#printf OUT "l 0 0 $w0z 0 0 $w1z\n";
        if ($dim == 3){
            printf OUT "l -$L2 -$L2 $w0z  $L2 -$L2 $w0z  \n";
            printf OUT "l  $L2 -$L2 $w0z  $L2  $L2 $w0z \n";
            printf OUT "l  $L2  $L2 $w0z -$L2  $L2 $w0z \n";
            printf OUT "l -$L2  $L2 $w0z -$L2 -$L2 $w0z  \n";
            printf OUT "l  0  $L2 $w0z  0 -$L2 $w0z  \n";
        }

        if ($dim == 3){
            printf OUT "l -$L2 -$L2 $w1z  $L2 -$L2 $w1z \n";
            printf OUT "l  $L2 -$L2 $w1z  $L2  $L2 $w1z \n";
            printf OUT "l  $L2  $L2 $w1z -$L2  $L2 $w1z \n";
			    printf OUT "l -$L2  $L2 $w1z -$L2 -$L2 $w1z \n";
            printf OUT "l   0   $L2 $w1z   0  -$L2 $w1z \n";
        }
		$lx2 = $lx/2;

		printf OUT "l  $lx2  0 $w0z -$lx2 0  $w0z \n";
		printf OUT "l  $lx2  0 $w1z -$lx2 0  $w1z \n";
		printf OUT "l  $lx2  0 $w0z  $lx2 0  $w1z \n";
		printf OUT "l -$lx2  0 $w0z -$lx2 0  $w1z \n";

		#		$rg = 1.3851 * ${vf}**(-1.97056);
		#printf OUT "r $rg\n";
		#printf OUT "@ 10\n";
		#printf OUT "c 0 0.01 0\n";
    }
}

#  RZ RY RX ==>  RY RZ RX
#	$rot1 =  atan2( 2*($q2*$q3 - $q0*$q1), $q0*$q0 - $q1*$q1 + $q2*$q2 - $q3*$q3);
#       $rot2 = -asin( 2*($q1*$q2 + $q0*$q3 ));
#      $rot3 = -atan2( 2*($q1*$q3 - $q0*$q2), $q0*$q0 + $q1*$q1 - $q2*$q2 - $q3*$q3);
#	$r1 =  -(180/3.1415926)*$rot1;
#	$r2 =  -(180/3.1415926)*$rot2;
#	$r3 =  -(180/3.1415926)*$rot3;
#	printf OUT "S(<$x, $z, $y>, <$r1, $r2, $r3>)\n";

#	$xx = $x;
#	$yy = $z;
#	$zz = $y;
#	$rxd = -(180/3.1415926)*$r1;
#	$ryd = -(180/3.1415926)*$r2;
#	$rzd = -(180/3.1415926)*$r3;
#	printf OUT "S (<$xx, $yy, $zz>, <$rxd, $ryd, $rzd>)\n";

sub printHead
{
    printf OUT 
        "#include \"colors.inc\"
#include \"shapes.inc\"
#include \"skies.inc\"
#include \"glass.inc\"
#include \"woods.inc\"
#include \"stones.inc\"
#include \"metals.inc\"
";

	printf OUT "
#macro S(p, rot)
sphere{ 0, 1
  pigment{checker White*1.2, Yellow scale 1}
  rotate rot
  translate p
 }
#end
";

    printf OUT 
        "
camera {
  location <0, 0, -$ARGV[1]>
  sky      <0, 1, 0>
  look_at  <0, 0, 0>
  angle 60
}
";

    printf OUT 
        "
light_source{
  <0, 100 , -100>
  color rgb <1, 1, 1>
  parallel
  point_at <0, 0, 0>
}
light_source{
  <0, 10 , -10>
  color rgb <1, 1, 1>
  parallel
  point_at <0, 0, 0>
}

light_source{
  <0, 0 , -100>
  color rgb <1, 1, 1>
  parallel
  point_at <0, 0, 0>
}

background {color rgb <0.02, 0.08, 0.15>}
";
}

sub yaplotColor {
	printf OUT "\@0 0 0 0 \n";
	#printf OUT "\@1 50 100 205 \n";
	printf OUT "\@1 25 50 102 \n";
	#printf OUT "\@1 255 255 255 \n";
	printf OUT "\@2 200 200 200 \n";
	printf OUT "\@3 50 150 255 \n";
	printf OUT "\@4 50 200 50 \n";
	printf OUT "\@5 255 100 100 \n";
	printf OUT "\@6 50 200 50 \n";
	printf OUT "\@7 255 255 0 \n";
	printf OUT "\@8 255 255 255\n";
	printf OUT "\@9 150 150 150\n";
	#printf OUT "\@8 224 143 0 \n";
	#printf OUT "\@9 67 163 230 \n";
	#printf OUT "\@8 253 105 6 \n";
	#printf OUT "\@9 109 109 109 \n";
	printf OUT "\@10 250 250 250 \n";
	printf OUT "\@11 240 240 240 \n";
	printf OUT "\@12 230 230 230 \n";
	printf OUT "\@13 220 220 220 \n";
	printf OUT "\@14 210 210 210 \n";
	printf OUT "\@15 200 200 200 \n";
	printf OUT "\@16 190 190 190 \n";
	printf OUT "\@17 180 180 180 \n";
	printf OUT "\@18 170 170 170 \n";
	printf OUT "\@19 160 160 160 \n";
	printf OUT "\@20 150 150 150 \n";
	printf OUT "\@21 140 140 140 \n";
	printf OUT "\@22 130 130 130 \n";
	printf OUT "\@23 120 120 120 \n";
	printf OUT "\@24 110 110 110 \n";
	printf OUT "\@25 100 100 100 \n";
	printf OUT "\@26 90 90 90 \n";
	printf OUT "\@27 80 90 90 \n";
	printf OUT "\@28 70 70 70\n";
	printf OUT "\@29 60 60 60 \n";
	printf OUT "\@30 50 50 50 \n";
	printf OUT "\@31 40 40 40 \n";
	printf OUT "\@32 30 30 30 \n";
	printf OUT "\@33 20 20 20 \n";
	printf OUT "\@34 10 10 10 \n";
	printf OUT "\@35 0 0 0 \n";
}


#object {
#  Plane_XZ
# translate -200*y
#  texture {T_Chrome_5C}
#}


#box{ <-50,-1,$wz0 ><50,1, $wz00 >
#material { texture { pigment { color Clear } finish { F_Glass1 } } interior { I_Glass1 fade_color Col_Emerald_03 } }
#}

#box{ <-50,-1,$wz1 ><50,1, $wz11 >
#material { texture { pigment { color Clear } finish { F_Glass1 } } interior { I_Glass1 fade_color Col_Emerald_03 } }
#}


#polygon{
#       5,
#       <-50,-1,$wz0 >, <50,-1,$wz0 >, <50,1,$wz0 >, <-50,1,$wz0 >, <-50,-1,$wz0 >
#       pigment{color rgb<0.8,0.8,0.4,0.5>*2}
#     }
#polygon{
#       5,
#       <-50,-1,$wz1 >, <50,-1,$wz1 >, <50,1,$wz1 >, <-50,1,$wz1 >, <-50,-1,$wz1>
#       pigment{color rgb<0.8,0.8,0.4,0.5>*2}
#     }


#material { texture { pigment { color Clear } finish { F_Glass1 } } interior { I_Glass1 fade_color Col_Emerald_03 } }
#pigment {color LightBlue}
# texture {T_Brass_5B}

#  texture{
#    pigment{color White}
#    finish{
#      phong 0.4
#      reflection 0.2
#    }
#  }

#printf OUT "
#camera {
#  location <0, -80, -80>
#  sky      <0, -1, 0>
#  look_at  <0, 0, -10>
#  angle 65
#}
#";



#printf OUT 
#"
##macro S(A, B)
#cylinder{ A, B, 0.05
#  pigment{color Red}
#}
##end
#";
