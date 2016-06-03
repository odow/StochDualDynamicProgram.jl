# The following a modified version of that found at
#
# Humanize.jl    https://github.com/IainNZ/Humanize.jl
# as it was at commit be2c55008b501e17ed13c0a9aa791d40214385ea
#
# Copyright (c) 2016 Oscar Dowson, Iain Dunning, Julian Gehring
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# O.D. 2016 renamed
const suffix = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]
# O.D. fix base
const base   = 1000.0

# O.D. 2016 remaned. drop style optoin
function humanize(value::Number, format="5.1f")
    # O.D. fix suffix
    # O.D. fix base
    bytes   = abs(float(value)) # O.D. abs value
    format  = "%$(format)%s"    # O.D. add % char to beginning
    fmt_str = @eval (v,s)->@sprintf($format,v,s)
    unit    = base
    s       = suffix[1]
    for (i,s) in enumerate(suffix)
        unit = base ^ (i)
        bytes < unit && break
    end
    # O.D. add sign
    return fmt_str(sign(value)*base * bytes / unit, s)::ASCIIString
end
