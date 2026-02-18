from django.core.mail import send_mail
from drugstone.settings import settings

default_sender = settings.EMAIL_ADDRESS


def bugreport(title, body, cc=None):
    if cc is None:
        send(title, body)
    else:
        send(title=title, body=body, recipient=['contact@drugst.one', cc])


def send(title, body, sender=default_sender, recipient=['contact@drugst.one'], fail_silently=False):
    send_mail(title, body, sender, recipient, fail_silently=fail_silently)
