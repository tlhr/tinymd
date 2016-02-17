let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Documents/code/test/rust/md/src
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +9 ~/Documents/code/test/rust/md/src/main.rs
badd +155 ~/Documents/code/theo/homework2/newmd.py
badd +22 ~/Documents/code/test/rust/md/src/linalg/linalg.rs
badd +10 ~/Documents/code/test/rust/md/Cargo.toml
badd +8 ~/Documents/code/test/rust/md/src/linalg/Cargo.toml
badd +1 ~/Documents/code/test/rust/pi/src/main.rs
badd +3 ~/Documents/code/pharmtools/rust-pepper/src/main.rs
argglobal
silent! argdel *
argadd main.rs
edit ~/Documents/code/test/rust/md/src/main.rs
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 89 + 89) / 178)
exe 'vert 2resize ' . ((&columns * 88 + 89) / 178)
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 40 - ((34 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
40
normal! 024|
wincmd w
argglobal
edit ~/Documents/code/test/rust/md/src/linalg/linalg.rs
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 181 - ((39 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
181
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 89 + 89) / 178)
exe 'vert 2resize ' . ((&columns * 88 + 89) / 178)
tabnext 1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
