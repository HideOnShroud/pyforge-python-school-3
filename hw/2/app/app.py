def lambda_handler(event, context):
    names = event.get('names', ['World'])
    greetings = [f'Hello, {name}!' for name in names]
    return {
        'statusCode': 200,
        'body': greetings
    }
