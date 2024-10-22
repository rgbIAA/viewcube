;Single Cube Califa Binaural Explorer
;Adrián García Riber
;2023

<CsoundSynthesizer>
<CsOptions>
-odac
</CsOptions>
<CsInstruments>
sr = 48000
ksmps = 32
nchnls = 2
0dbfs = 1

gks init 0

gkamp0 init 0
gkamp1 init 0
gkamp2 init 0
gkamp3 init 0
gkamp4 init 0
gkamp5 init 0

gkf0 init 0
gkf1 init 0
gkf2 init 0
gkf3 init 0
gkf4 init 0
gkf5 init 0

gkA0 init 1
gkA1 init 1
gkA2 init 1
gkA3 init 1
gkA4 init 1
gkA5 init 1

gkaz init 0
gkdist init 0
gkflux init 0

gkFader init 0.5
gkLpFrec init 20000
gkHpFrec init 20

chn_S "hrtf_L", 1
chn_S "hrtf_R", 1



instr 1

giosc_amps OSCinit 9970
giosc_freqs OSCinit 9971
giosc_coords OSCinit 9972
giosc_flux OSCinit 9973

gS_HRTF_left chnget "hrtf_L"
gS_HRTF_right chnget "hrtf_R"

iosc_s OSCinit 9999

kans_amps OSClisten giosc_amps, "amp0/amp1/amp2/amp3/amp4/amp5", "ffffff", gkamp0 ,gkamp1, gkamp2, gkamp3, gkamp4,gkamp5
;printk2 gkamp0
;printk2 gkamp1
;printk2 gkamp2
;printk2 gkamp3
;printk2 gkamp4
;printk2 gkamp5

kans_freqs OSClisten giosc_freqs, "lat0/lat1/lat2/lat3/lat4/lat5", "ffffff", gkf0 ,gkf1, gkf2, gkf3, gkf4,gkf5
;printk2 gkf0
;printk2 gkf1
;printk2 gkf2
;printk2 gkf3
;printk2 gkf4
;printk2 gkf5
 
kans_coords OSClisten giosc_coords, "az/dist", "ff", gkaz ,gkdist
;printk2 gkaz
;printk2 gkdist

kans_coords OSClisten giosc_flux, "flux", "f", gkflux
;printk2 gkflux


;Additive Synth

kAmp = gkflux

a0  oscil kAmp*gkamp0, gkf0, -1, 0
if gkf0<0 then
    a0 = -a0
endif
a1  oscil kAmp*gkamp1, gkf1, -1, 0
if gkf1<0 then
    a1 = -a1
endif
a2  oscil kAmp*gkamp2, gkf2, -1, 0
if gkf2<0 then
    a2 = -a2
endif
a3	oscil kAmp*gkamp3, gkf3, -1, 0
if gkf3<0 then
    a3 = -a3
endif
a4	oscil kAmp*gkamp4, gkf4, -1, 0
if gkf4<0 then
    a4 = -a4
endif
a5	oscil kAmp*gkamp5, gkf5, -1, 0
if gkf5<0 then
    a5 = -a5
endif

aOut = (a0+a1+a2+a3+a4+a5)/4

adel linseg 0, .5, 0.2, .5, 0	;max delay time =20ms	

aDly comb aOut, 30, .5
aDlyFilt moogladder aDly, 700, .1
aDlyFX phaser1 aDlyFilt, 10, 1, .5
aflg flanger aDlyFX, adel, .5, 10

aHp butterhp aOut-aflg/8, gkHpFrec
aFilt moogladder aHp, gkLpFrec, 0.1

gasendL= aFilt*gkdist	
gasendR= aFilt*gkdist

gaRevLf, gaRevRf		reverbsc	gasendL,gasendR,.9,10000

aLeft, aRight  hrtfmove   (aFilt-gaRevLf)*gkFader*gkflux, gkaz, 0, gS_HRTF_left, gS_HRTF_right, 4, 9.0, 48000

               outs        aLeft, aRight                 ; ---------------------------binaural outputs


endin


</CsInstruments>
<CsScore>
f 1 0 1024 10 1
i 1 0 605000
e
</CsScore>
</CsoundSynthesizer>
<bsbPanel>
 <label>Widgets</label>
 <objectName/>
 <x>100</x>
 <y>100</y>
 <width>320</width>
 <height>240</height>
 <visible>true</visible>
 <uuid/>
 <bgcolor mode="nobackground">
  <r>255</r>
  <g>255</g>
  <b>255</b>
 </bgcolor>
</bsbPanel>
<bsbPresets>
</bsbPresets>
