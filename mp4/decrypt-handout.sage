import struct
import re, time
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from base64 import b64decode
from math import ceil, floor
from pbp import int_to_mpi, parse_mpi, unpad


### Message formatting ###
encrypt_header = b'-----BEGIN PRETTY BAD ENCRYPTED MESSAGE-----\n'
encrypt_footer = b'-----END PRETTY BAD ENCRYPTED MESSAGE-----\n'
debug = 0

# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii,jj] == 0 else 'X'
            a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print (a)

# function derived from https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/coppersmith.sage
def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Coppersmith revisited by Howgrave-Graham
    
    finds a solution if:
    * b|modulus, b >= modulus^beta , 0 < beta <= 1
    * |x| < XX
    """
    #
    # init
    #
    dd = pol.degree()
    nn = dd * mm + tt

    #
    # checks
    #
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    #
    # calculate bounds and display them
    #
    """
    * we want to find g(x) such that ||g(xX)|| <= b^m / sqrt(n)
    * we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)
    * we will use that vector as a coefficient vector for our g(x)
    
    * so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
    
    so we can obtain ||v|| < N^(beta*m) / sqrt(n) <= b^m / sqrt(n)
    (it's important to use N because we might not know b)
    """
    if debug:
        # t optimized?
        print ("\n# Optimized t?\n")
        print ("we want X^(n-1) < N^(beta*m) so that each vector is helpful")
        cond1 = RR(XX^(nn-1))
        print ("* X^(n-1) = ", cond1)
        cond2 = pow(modulus, beta*mm)
        print ("* N^(beta*m) = ", cond2)
        print ("* X^(n-1) < N^(beta*m) \n-> GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) \n-> NOT GOOD")
        
        # bound for X
        print ("\n# X bound respected?\n")
        print ("we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M")
        print ("* X =", XX)
        cond2 = RR(modulus^(((2*beta*mm)/(nn-1)) - ((dd*mm*(mm+1))/(nn*(nn-1)))) / 2)
        print ("* M =", cond2)
        print ("* X <= M \n-> GOOD" if XX <= cond2 else "* X > M \n-> NOT GOOD")

        # solution possible?
        print ("\n# Solutions possible?\n")
        detL = RR(modulus^(dd * mm * (mm + 1) / 2) * XX^(nn * (nn - 1) / 2))
        print ("we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)")
        cond1 = RR(2^((nn - 1)/4) * detL^(1/nn))
        print ("* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1)
        cond2 = RR(modulus^(beta*mm) / sqrt(nn))
        print ("* N^(beta*m) / sqrt(n) = ", cond2)
        print ("* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) \n-> SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) \n-> NO SOLUTIONS MIGHT BE FOUND (but we never know)")

        # warning about X
        print ("\n# Note that no solutions will be found _for sure_ if you don't respect:\n* |root| < X \n* b >= modulus^beta\n")
    
    #
    # Coppersmith revisited algo for univariate
    #

    # change ring of pol and x
 
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    # display basis matrix
    # if debug:
    matrix_overview(BB, modulus^mm)

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()
    print ("potential roots:", potential_roots)

    # test roots
    roots = []
    for root in potential_roots:
        if root[0].is_integer():
            result = polZ(ZZ(root[0]))
            if gcd(modulus, result) >= modulus^beta:
                roots.append(ZZ(root[0]))

    return roots

def decrypt():
    # Read the file
    s = open('mp4-ece498ac-fall2019.enc.asc','rb').read()
    # s = "-----BEGIN PRETTY BAD ENCRYPTED MESSAGE-----\nAAEAAPTMPKyqs/WrCXWxUPM1uw2uWmUeIlRmNw+0YWtxLIk6xAm1UBiFwPMC+RZnZbUiPoMwcpk+MADyU5CWwMhKtG7CDB2JMqB5dJ2fJ3dMMU+5taTwa466dD+h4V3uwR/GB/ihxVsg2JVgT0/aiR4XXDDB1HI8g5YLeDMbzlg1AX7uT36YcU2MZEajymmjR94hpHg6qLNtZtYNDZE87x/gJWc0MEDd4td62xRdyImT885kW7jR+go8mmBXg2rLqOtrQXlGImZv2E2l72QoPIOQ94EvCF5A/+PBJqzZqruvuAR+CUuqhgSohHGP7/c8FFc+xwiS6cJAfnWjeCQnpj8wGgCsOfufkQ85ymmv9Gk+zWaE94d/y1nJJTxiff/gpRh/og==\n-----END PRETTY BAD ENCRYPTED MESSAGE-----\n"
    # s = s.encode()
    # Remove header
    data = re.search(encrypt_header+b"(.*)"+encrypt_footer,s,flags=re.DOTALL).group(1)
    data = b64decode(data)

    # Split into rsa encryption c and symmetric encryption (iv, cipher)
    c, index = parse_mpi(data,0)
    iv = data[index:index+AES.block_size]
    ciphertext = data[index+AES.block_size:]

    # Import the public key for coppersmith attack
    pubkey = RSA.importKey(open('key.pub').read())

    # Problem to equation (default)
    P.<x> = PolynomialRing(Zmod(pubkey.n))

    # multiply x with mask in order to account for the padding scheme k2 || k2 || k2 || k2
    mask = 1 + (1 << 256) + (1 << 256*2) + (1 << 256*3)
    poly = (x*mask)^pubkey.e - c
    # make the polynomial monic by dividing the coefficient of the highest order term
    poly /= mask^pubkey.e
    # print(pol)
    dd = poly.degree()

    # Tweak those
    # luckily, q >= N^beta
    beta = 0.5                              # b = 0.5 we know the higher bits of N apparently
    epsilon = beta / 7                      # <= beta / 7
    mm = ceil(beta**2 / (dd * epsilon))     # optimized value
    tt = floor(dd * mm * ((1/beta) - 1))    # optimized value
    XX = ceil(pubkey.n**((beta**2/dd) - epsilon)) + 1000000000000000000000000000000000 # optimized value

    # Coppersmith
    start_time = time.time()
    roots = coppersmith_howgrave_univariate(poly, pubkey.n, beta, mm, tt, XX)
    # output
    print ("\n# Solutions")
    print ("we found:", str(roots))
    print("in: %s seconds " % (time.time() - start_time))

    # length = 16,24 not enough to contain the integer size
    cipher = AES.new(int(roots[0]).to_bytes(length=32, byteorder='big')[:16], AES.MODE_CBC, iv)
    plaintext = unpad(cipher.decrypt(ciphertext))
    # print(plaintext)
    f = open("mp4-decode.pdf",'wb')
    f.write(plaintext)
    f.close()


if __name__ == "__main__":
    decrypt()
