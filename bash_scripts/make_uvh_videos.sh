python plotting.py $1

cd $1/plots/

#ffmpeg -f image2 -framerate 24 -i u_%04d.png u.mp4
#ffmpeg -f image2 -framerate 24 -i v_%04d.png v.mp4
ffmpeg -f image2 -framerate 24 -i h_%04d.png h.mp4

rm *_*
