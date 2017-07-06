#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"

#include <cstdio>
#include <cassert>
#include <unistd.h>

#include <stdexcept>
#include <fcntl.h>

bool CloseCoutSentry::open_ = true;
int  CloseCoutSentry::fdOut_ = 0;
int  CloseCoutSentry::fdErr_ = 0;
int  CloseCoutSentry::fdTmp_ = 0;
int  CloseCoutSentry::fdOutDup_ = 0;
FILE * CloseCoutSentry::trueStdOut_ = 0;
CloseCoutSentry *CloseCoutSentry::owner_ = 0;


CloseCoutSentry::CloseCoutSentry(bool silent) :
    silent_(silent), stdOutIsMine_(false)
{
    if (silent_) {
        if (open_) {
            open_ = false;
            if (fdOut_ == 0 && fdErr_ == 0) {
                fdOut_ = dup(1);
                fdErr_ = dup(2);
            }
            fdTmp_ = open( "/dev/null", O_RDWR );
            dup2(fdTmp_, 1);
            dup2(fdTmp_, 2);
            assert(owner_ == 0);
            owner_ = this;
        } else {
            silent_ = false; 
        }
    }
}

CloseCoutSentry::~CloseCoutSentry() 
{
    clear();
}

void CloseCoutSentry::clear() 
{
    if (stdOutIsMine_) { 
        assert(this == owner_);
        fclose(trueStdOut_); trueStdOut_ = 0; stdOutIsMine_ = false;
    }
    if (silent_) {
        reallyClear();
        silent_ = false;
    }
}

void CloseCoutSentry::reallyClear() 
{
    if (fdOut_ != fdErr_) {
        dup2( fdOut_, 1 );
        dup2( fdErr_, 2 );
        open_   = true;
        owner_ = 0;
    }
}

void CloseCoutSentry::breakFree() 
{
    reallyClear();
}

FILE *CloseCoutSentry::trueStdOutGlobal()
{
    if (!owner_) return stdout;
    return owner_->trueStdOut();
}

FILE *CloseCoutSentry::trueStdOut() 
{
    if (open_) return stdout;
    if (trueStdOut_) return trueStdOut_;
    if (owner_ != this && owner_ != 0) return owner_->trueStdOut();
    assert(owner_ == this);
    stdOutIsMine_ = true;
    fdOutDup_ = dup( fdOut_ ); // When clear() calls fclose(trueStdOut_), this makes sure fdOutDup_ gets closed instead of fdOut_
    trueStdOut_ = fdopen( fdOutDup_, "w" );
    return trueStdOut_;
}
