import smtplib
from email.mime.text import MIMEText

def send(subject,body,recipient='ashbaker@sas.upenn.edu',attachment=None):
    # contact here
    recipient  = recipient
    # login credentials ...should probably worry about this
    username = 'camal.telescope'
    password = 'C@m@l20!6'
    # Prep email
    msg = MIMEText(body)
    msg['From']    = username
    msg['To']      = recipient
    msg['Subject'] = subject
    # send email
    server = smtplib.SMTP('smtp.gmail.com')
    server.starttls()
    server.login(username,password)
    server.sendmail(username, recipient, msg.as_string())
    server.quit()

if __name__ == "__main__":
    send('test subject','body of email')