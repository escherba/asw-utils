/*
 * =====================================================================================
 *
 *       Filename:  BamIterator.hpp
 *
 *    Description:  A wrapper around BAM reading functions provided in samtools
 *
 *        Version:  1.0
 *        Created:  04/08/2011 15:03:18
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

class BamIterator;
/*
 * =====================================================================================
 *        Class:  BamHeader
 *  Description:  A wrapper around Samtools bam_header_t type
 * =====================================================================================
 */
class BamHeader
{
private:
        bam_header_t *header;
        // disallow copying and asignment by declaring the copy constructor and the
        // assignment operator private
        BamHeader(const BamHeader&);
        BamHeader& operator=(const BamHeader&);
public:
/*
 *--------------------------------------------------------------------------------------
 *       Class:  BamHeader
 *      Method:  BamHeader :: BamHeader
 * Description:  Default constructor
 *--------------------------------------------------------------------------------------
 */
        BamHeader(bamFile fp) :
                header(bam_header_read(fp))
        {}

        BamHeader() :
                header(NULL)
        {}
/*
 *--------------------------------------------------------------------------------------
 *       Class:  BamHeader
 *      Method:  BamHeader :: ~BamHeader
 * Description:  Default destructor
 *--------------------------------------------------------------------------------------
 */
        ~BamHeader() {
                if (header) bam_header_destroy(header);
        }
        bam_header_t* operator->() { return header; }
        bam_header_t& operator*() { return *header; }

        friend class BamIterator;
};


class BamIterator : public std::iterator<std::input_iterator_tag, bam1_t>
{
public:
        // Copy Constuctor
        BamIterator(const BamIterator& r) : 
                _is_owner(false), 
                _is_fp_owner(false),
                _row_number(r._row_number),
                _fp(r._fp), 
                _bam(r._bam), 
                _read_return(r._read_return) 
        {}
        BamIterator() : 
                _is_owner(false), 
                _is_fp_owner(false),
                _row_number(0u),
                _fp(NULL), 
                _bam(NULL), 
                _read_return(0) 
        {}
        BamIterator(char *path, BamHeader& bh) :
                _is_owner(true),
                _is_fp_owner(true),
                _row_number(1u),
                _fp(bam_open(path, "r")),
                _bam(bam_init1()),
                _read_return(0)
        {
                if (_fp == NULL) {
                        bam_destroy1(_bam);
                        exit(EXIT_FAILURE);
                }
                bh.header = bam_header_read(_fp);
                _read_return = bam_read1(_fp, _bam);
                if (_read_return <= 0) {
                        _row_number = 0u;
                }
        }
        BamIterator(bamFile fp) : 
                _is_owner(true), 
                _is_fp_owner(false),
                _row_number(1u),
                _fp(fp), 
                _bam(bam_init1()), 
                _read_return(bam_read1(fp, _bam))
        {
                if (_read_return <= 0) {
                        _row_number = 0u;
                }
        }
        ~BamIterator() {
                if (_is_fp_owner) bam_close(_fp);
                if (_is_owner) bam_destroy1(_bam);
        }
        // --------------------------------
        BamIterator& operator=(const BamIterator& r) {
                // Treating self-asignment (this == &r) the same as assignment
                _is_owner = false;
                _is_fp_owner = r._is_fp_owner;
                _row_number = r._row_number;
                _fp = r._fp;
                _bam = r._bam;
                _read_return = r._read_return;
                return *this;
        }
        BamIterator& operator++() { // pre-increment
                _read_return = bam_read1(_fp, _bam);
                if (_read_return > 0) {
                        ++_row_number;
                } else {
                        std::cerr << _row_number << " rows read" << std::endl;
                        _row_number = 0u;
                }
                return *this;
        }
        BamIterator operator++(int) { // post-increment
                BamIterator tmp(*this);
                operator++();
                return tmp;
        }
        bool operator==(const BamIterator& r) {
                return _row_number == r._row_number;
        }
        bool operator!=(const BamIterator& r) {
                return _row_number != r._row_number;
        }
        bam1_t* operator->() { return _bam; }
        bam1_t& operator*() { return *_bam; }

protected:
        bool _is_owner; // whether to call bam_destroy1 upon object destruction
        bool _is_fp_owner;
        unsigned int _row_number;
        bamFile _fp;
        bam1_t* _bam;
        int _read_return;
};


