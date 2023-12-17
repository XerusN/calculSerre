set terminal gif animate delay 1

set cbrange[260:600]

set output '~/Documents/VSCode/hugeFile/serre/serre1.gif'

do for [i = 0:500] {
    plot '~/Documents/VSCode/hugeFile/serre/serreFortranOutput.txt' index i matrix w image
}