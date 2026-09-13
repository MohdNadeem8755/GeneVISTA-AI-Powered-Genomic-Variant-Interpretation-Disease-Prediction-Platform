"""Run locally; paste only the resulting salted hash into Render secrets."""
import getpass
import hashlib
import secrets
if __name__=='__main__':
    password=getpass.getpass('Workspace password (at least 12 characters): ')
    if len(password)<12: raise SystemExit('Use at least 12 characters.')
    if password!=getpass.getpass('Confirm password: '): raise SystemExit('Passwords do not match.')
    salt=secrets.token_bytes(16)
    print(salt.hex()+':'+hashlib.scrypt(password.encode(),salt=salt,n=16384,r=8,p=1).hex())
