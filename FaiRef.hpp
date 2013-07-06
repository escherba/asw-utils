/*
 * =====================================================================================
 *
 *       Filename:  FaiRef.hpp
 *
 *    Description:  A wrapper to load genome reference contigs from FAI-indexed FASTA
 *                  files.
 *
 *        Version:  1.0
 *        Created:  04/08/2011 15:06:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eugene Scherba (es), escherba@bu.edu
 *        Company:  Boston University
 *
 * Copyright (c) 2011 Eugene Scherba
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * =====================================================================================
 */

/*
 * =====================================================================================
 *        Class:  FaiRef
 *  Description:  
 * =====================================================================================
 */
class FaiRef
{
public:
        /* ====================  LIFECYCLE     ======================================= */
        FaiRef (const char *path, BamHeader& cbh) :
                fai(fai_load(path)),
                header(&cbh),
                tid(-2),
                ref(NULL),
                ref_len(0)
        {};                             /* constructor */
        ~FaiRef() {
                fai_destroy(fai);
                if (ref) free(ref);
        }
        /* ====================  ACCESSORS     ======================================= */
        /* ====================  MUTATORS      ======================================= */
        template <typename T>
        const char* loadContig(int new_tid, T offset = 0) throw(char*) {
                if (new_tid != tid) {
                        if (ref) free(ref);
                        tid = new_tid;
                        ref = fai_fetch(fai, (*header)->target_name[tid], &ref_len);
                }
                return getCoord<T>(offset);
        }
        template <typename T>
        const char* getCoord(T offset) throw(char*) {
                // Bounds the coordinate to a valid interval within reference sequence
                if (ref == NULL) {
                        throw((*header)->target_name[tid]);
                        //std::cerr << "ERROR: failed to find sequence " 
                        //          << (*header)->target_name[tid] 
                        //          << "  in the reference." 
                        //          << std::endl;
                }
                return ref + std::min<T>(ref_len, std::max<T>(0, offset));
        }
        /* ====================  OPERATORS     ======================================= */

protected:
        /* ====================  DATA MEMBERS  ======================================= */

private:
        faidx_t *fai;
        BamHeader *header;
        int tid;
        char *ref;
        int ref_len;

        FaiRef(const FaiRef&);
        FaiRef& operator=(const FaiRef&);
        /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class FaiRef  ----- */


